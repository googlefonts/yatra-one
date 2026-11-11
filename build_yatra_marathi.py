import os
import sys
import copy
from fontTools.ttLib import TTFont
from fontTools.ttLib.tables import otTables
from fontTools.ttLib.tables._g_l_y_f import Glyph, GlyphComponent
import uharfbuzz as hb


def add_headline_bar(font, glyph_name, headline_y, bar_thickness=45):
    glyf = font['glyf']
    if glyph_name not in glyf:
        return False
    g = glyf[glyph_name]
    if g.numberOfContours <= 0 or g.isComposite():
        return False

    xMin, yMin, xMax, yMax = g.xMin, g.yMin, g.xMax, g.yMax
    if (xMax - xMin) < 100:
        return False

    bar_top = yMax
    bar_bottom = yMax - bar_thickness

    bar_coords = [
        (xMin, bar_top),
        (xMax, bar_top),
        (xMax, bar_bottom),
        (xMin, bar_bottom),
    ]

    coords = list(g.coordinates)
    flags = list(g.flags)
    endPts = list(g.endPtsOfContours)

    new_start = len(coords)
    for pt in bar_coords:
        coords.append(pt)
        flags.append(1)
    endPts.append(new_start + len(bar_coords) - 1)

    g.coordinates = type(g.coordinates)(coords)
    g.flags = type(g.flags)(flags)
    g.endPtsOfContours = endPts
    g.numberOfContours = len(endPts)
    g.recalcBounds(glyf)
    return True


def apply_headline_fixes(font):
    glyf = font['glyf']
    sample = ['dvKA', 'dvGA', 'dvTA', 'dvDA', 'dvNA', 'dvPA', 'dvMA', 'dvSA']
    heights = [glyf[name].yMax for name in sample if name in glyf and glyf[name].numberOfContours > 0]
    headline_y = (sum(heights) / len(heights)) if heights else 750.0

    conjuncts = [
        name for name in font.getGlyphOrder()
        if name.startswith('dv') and '_' in name and not name.startswith('dvm')
    ]
    modified = 0
    for name in conjuncts:
        if add_headline_bar(font, name, headline_y, bar_thickness=45):
            modified += 1
    print(f"Headline bars added to {modified} conjunct glyphs.")


def extract_ligatures(lookup):
    rules = []
    if lookup.LookupType != 4:
        return rules
    for st in lookup.SubTable:
        if hasattr(st, 'ligatures'):
            for fg, ligs in st.ligatures.items():
                for l in ligs:
                    rules.append(([fg] + list(l.Component), l.LigGlyph))
    return rules


def build_lig_subtable(rules):
    st = otTables.LigatureSubst()
    st.Format = 1
    st.ligatures = {}
    for comps, result in rules:
        first = comps[0]
        lig = otTables.Ligature()
        lig.LigGlyph = result
        lig.Component = comps[1:]
        lig.CompCount = len(comps)
        if first not in st.ligatures:
            st.ligatures[first] = []
        st.ligatures[first].append(lig)

    for fg in st.ligatures:
        st.ligatures[fg].sort(key=lambda l: -l.CompCount)
    return st


def apply_conjunct_and_marathi_rules(font):
    gsub = font['GSUB'].table
    virama = 'dvVirama'

    half_map = {}
    for li in [11, 12]:
        for comps, res in extract_ligatures(gsub.LookupList.Lookup[li]):
            if len(comps) == 2 and comps[1] == virama:
                half_map[res] = comps[0]

    pres_rules = extract_ligatures(gsub.LookupList.Lookup[14])
    rkrf_rules = extract_ligatures(gsub.LookupList.Lookup[8])
    akhn_orig = extract_ligatures(gsub.LookupList.Lookup[5])

    three_glyph = []
    for comps, result in pres_rules:
        if len(comps) == 2:
            hf, base = comps
            if hf in half_map:
                three_glyph.append(([half_map[hf], virama, base], result))
        elif len(comps) == 3:
            h1, h2, base = comps
            if h1 in half_map and h2 in half_map:
                three_glyph.append(([half_map[h1], virama, half_map[h2], virama, base], result))
            elif h1 in half_map:
                three_glyph.append(([half_map[h1], virama, h2, base], result))
        elif len(comps) == 4:
            h1, h2, h3, base = comps
            if h1 in half_map and h2 in half_map and h3 in half_map:
                three_glyph.append(([half_map[h1], virama, half_map[h2], virama, half_map[h3], virama, base], result))

    existing = set(tuple(r[0]) for r in three_glyph)
    for comps, result in rkrf_rules:
        if tuple(comps) not in existing:
            three_glyph.append((comps, result))
            existing.add(tuple(comps))

    for comps, result in akhn_orig:
        if tuple(comps) not in existing:
            three_glyph.append((comps, result))
            existing.add(tuple(comps))

    mar_variants = []
    for comps, result in three_glyph:
        if 'dvLA' in comps or 'dvSHA' in comps:
            new_comps = [
                'dvLA.mar' if c == 'dvLA' else ('dvSHA.mar' if c == 'dvSHA' else c)
                for c in comps
            ]
            if tuple(new_comps) not in existing:
                mar_variants.append((new_comps, result))
                existing.add(tuple(new_comps))

    three_glyph.extend(mar_variants)

    half_rules_11 = extract_ligatures(gsub.LookupList.Lookup[11])
    half_rules_11.append((['dvLA.mar', virama], 'dvL'))
    half_rules_11.append((['dvSHA.mar', virama], 'dvSH'))

    merged_rules = three_glyph + half_rules_11

    new_lookup = otTables.Lookup()
    new_lookup.LookupType = 4
    new_lookup.LookupFlag = 0
    new_lookup.SubTable = [build_lig_subtable(merged_rules)]

    new_half_idx = len(gsub.LookupList.Lookup)
    gsub.LookupList.Lookup.append(new_lookup)
    gsub.LookupList.LookupCount = len(gsub.LookupList.Lookup)

    for fr in gsub.FeatureList.FeatureRecord:
        if fr.FeatureTag == 'half':
            old = list(fr.Feature.LookupListIndex)
            new_list = [new_half_idx if x == 11 else x for x in old]
            fr.Feature.LookupListIndex = new_list
            fr.Feature.LookupCount = len(new_list)

    akhn_lookup = otTables.Lookup()
    akhn_lookup.LookupType = 4
    akhn_lookup.LookupFlag = 0
    akhn_lookup.SubTable = [build_lig_subtable(three_glyph)]
    akhn_idx = len(gsub.LookupList.Lookup)
    gsub.LookupList.Lookup.append(akhn_lookup)
    gsub.LookupList.LookupCount = len(gsub.LookupList.Lookup)

    for fr in gsub.FeatureList.FeatureRecord:
        if fr.FeatureTag == 'akhn':
            fr.Feature.LookupListIndex = [akhn_idx]
            fr.Feature.LookupCount = 1

    cjct_lookup = otTables.Lookup()
    cjct_lookup.LookupType = 4
    cjct_lookup.LookupFlag = 0
    cjct_lookup.SubTable = [build_lig_subtable(three_glyph)]
    cjct_idx = len(gsub.LookupList.Lookup)
    gsub.LookupList.Lookup.append(cjct_lookup)
    gsub.LookupList.LookupCount = len(gsub.LookupList.Lookup)

    cjct_rec = otTables.FeatureRecord()
    cjct_rec.FeatureTag = 'cjct'
    cjct_rec.Feature = otTables.Feature()
    cjct_rec.Feature.FeatureParams = None
    cjct_rec.Feature.LookupListIndex = [cjct_idx]
    cjct_rec.Feature.LookupCount = 1
    gsub.FeatureList.FeatureRecord.append(cjct_rec)
    gsub.FeatureList.FeatureCount = len(gsub.FeatureList.FeatureRecord)
    cjct_feat = gsub.FeatureList.FeatureCount - 1

    for sr in gsub.ScriptList.ScriptRecord:
        if sr.ScriptTag in ('DFLT', 'dev2', 'deva'):
            if sr.Script.DefaultLangSys:
                if cjct_feat not in sr.Script.DefaultLangSys.FeatureIndex:
                    sr.Script.DefaultLangSys.FeatureIndex.append(cjct_feat)
                    sr.Script.DefaultLangSys.FeatureCount = len(sr.Script.DefaultLangSys.FeatureIndex)
                dflt_features = list(sr.Script.DefaultLangSys.FeatureIndex)
            else:
                dflt_features = []

            for lr in sr.Script.LangSysRecord:
                combined = list(lr.LangSys.FeatureIndex)
                for fi in dflt_features:
                    if fi not in combined:
                        combined.append(fi)
                lr.LangSys.FeatureIndex = combined
                lr.LangSys.FeatureCount = len(combined)
    print("Marathi LangSys MAR connected with all Devanagari features.")


def apply_aa_bindi_fixes(font):
    glyf = font['glyf']
    hmtx = font['hmtx']
    gsub = font['GSUB'].table

    comp_aa = GlyphComponent()
    comp_aa.glyphName = 'dvmAA'
    comp_aa.x = 0
    comp_aa.y = 0
    comp_aa.flags = 0x0204

    comp_anu = GlyphComponent()
    comp_anu.glyphName = 'dvAnusvara'
    comp_anu.x = 55
    comp_anu.y = 0
    comp_anu.flags = 0

    g_lig = Glyph()
    g_lig.numberOfContours = -1
    g_lig.components = [comp_aa, comp_anu]
    glyf['dvmAA_Anusvara'] = g_lig
    hmtx['dvmAA_Anusvara'] = (269, -5)

    comp_can = GlyphComponent()
    comp_can.glyphName = 'dvCandrabindu'
    comp_can.x = 320
    comp_can.y = 0
    comp_can.flags = 0

    g_lig_can = Glyph()
    g_lig_can.numberOfContours = -1
    g_lig_can.components = [comp_aa, comp_can]
    glyf['dvmAA_Candrabindu'] = g_lig_can
    hmtx['dvmAA_Candrabindu'] = (269, -5)

    glyph_order = font.getGlyphOrder()
    for gname in ['dvmAA_Anusvara', 'dvmAA_Candrabindu']:
        if gname not in glyph_order:
            glyph_order.append(gname)
    font.setGlyphOrder(glyph_order)

    lookup17 = gsub.LookupList.Lookup[17]
    subtable = lookup17.SubTable[0]

    lig_anu = otTables.Ligature()
    lig_anu.Component = ['dvAnusvara']
    lig_anu.LigGlyph = 'dvmAA_Anusvara'

    lig_can = otTables.Ligature()
    lig_can.Component = ['dvCandrabindu']
    lig_can.LigGlyph = 'dvmAA_Candrabindu'

    if 'dvmAA' not in subtable.ligatures:
        subtable.ligatures['dvmAA'] = []

    subtable.ligatures['dvmAA'].append(lig_anu)
    subtable.ligatures['dvmAA'].append(lig_can)
    print("Issue #5 (Aa + bindi centering) resolved.")


def apply_marathi_matra_and_anusvara_fixes(font):
    gsub = font['GSUB'].table
    gpos = font['GPOS'].table

    # 1. GSUB Lookup 15 (pres: dvmI -> dvmI.a04 for dvSHA.mar and dvLA.mar)
    l15 = gsub.LookupList.Lookup[15]
    for st in l15.SubTable:
        if hasattr(st, 'LookAheadCoverage') and st.LookAheadCoverage:
            cov = st.LookAheadCoverage[0].glyphs
            if 'dvSHA' in cov and 'dvSHA.mar' not in cov:
                cov.append('dvSHA.mar')
            if 'dvLA' in cov and 'dvLA.mar' not in cov:
                cov.append('dvLA.mar')

    # 2. GSUB Lookup 18 (abvs: dvAnusvara -> dvAnusvara.amI after dvmI.aXX + dvSHA.mar/dvLA.mar)
    l18 = gsub.LookupList.Lookup[18]
    for st in l18.SubTable:
        if hasattr(st, 'BacktrackCoverage') and st.BacktrackCoverage:
            for b_cov in st.BacktrackCoverage:
                cov = b_cov.glyphs
                if 'dvSHA' in cov and 'dvSHA.mar' not in cov:
                    cov.append('dvSHA.mar')
                if 'dvLA' in cov and 'dvLA.mar' not in cov:
                    cov.append('dvLA.mar')

    # 3. GSUB Lookup 19 (blws: dvSHA.mar + dvmvR -> dvSHA_mvR)
    l19 = gsub.LookupList.Lookup[19]
    for st in l19.SubTable:
        if hasattr(st, 'ligatures'):
            if 'dvSHA' in st.ligatures and 'dvSHA.mar' not in st.ligatures:
                st.ligatures['dvSHA.mar'] = []
                for lig in st.ligatures['dvSHA']:
                    new_lig = copy.deepcopy(lig)
                    st.ligatures['dvSHA.mar'].append(new_lig)

    # 4. GPOS Lookup 2 (abvm: anchor for dvAnusvara.amI on dvSHA.mar and dvLA.mar)
    l2 = gpos.LookupList.Lookup[2]
    for st in l2.SubTable:
        base_cov = st.BaseCoverage.glyphs
        base_array = st.BaseArray.BaseRecord
        for g_src, g_target in [('dvSHA', 'dvSHA.mar'), ('dvLA', 'dvLA.mar')]:
            if g_src in base_cov and g_target not in base_cov:
                src_idx = base_cov.index(g_src)
                src_rec = base_array[src_idx]
                base_cov.append(g_target)
                base_array.append(copy.deepcopy(src_rec))
    print("Injected Marathi matra & anusvara rules for dvSHA.mar and dvLA.mar.")


def update_yatra_metadata(font):
    name_table = font['name']
    name_table.setName("Version 1.002; 2026; Marathi OpenType shaping bug fixes", 5, 3, 1, 0x409)
    name_table.setName("Version 1.002; 2026; Marathi OpenType shaping bug fixes", 5, 1, 0, 0)
    font['head'].fontRevision = 1.002


def build_file(file_path):
    print(f"\nProcessing {file_path}...")
    font = TTFont(file_path)
    apply_headline_fixes(font)
    apply_conjunct_and_marathi_rules(font)
    apply_aa_bindi_fixes(font)
    apply_marathi_matra_and_anusvara_fixes(font)
    update_yatra_metadata(font)
    font.save(file_path)
    print(f"Saved {file_path}")


if __name__ == '__main__':
    base_dir = r"d:\web\projects\yatra-one"
    for fname in ["Yatra-Regular.ttf", "Yatra-Regular-Hinted.ttf"]:
        p = os.path.join(base_dir, fname)
        if os.path.exists(p):
            build_file(p)
