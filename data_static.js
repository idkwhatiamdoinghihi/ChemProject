/*
 ==========================================================================
 ★ Visual Bond Length Standards v16.2
 ==========================================================================
 Reference: anchored to 1,2-dichloropropane (C-C ~ 70, C-H ~ 50, C-Cl ~ 75)
 
 [1] Base radius contribution
 --------------------------------------------------------------------------
  - H (Hydrogen) .............. 15  (smallest, keeps structures compact)
  - Row 2 (C, N, O, F) ........ 35  (baseline)
  - Row 3 (Si, P, S, Cl) ...... 40  (slightly larger)
  - Row 4 (Br) ................ 45
  - Row 5 (I, Xe) ............. 50  (largest)

 [2] Bond order multipliers
 --------------------------------------------------------------------------
  - Single bond ................ x 1.00
  - Double bond ................ x 0.90
  - Triple bond ................ x 0.85

 [3] Example calculations (common bond lengths)
 --------------------------------------------------------------------------
  Type      Calc (R1 + R2) * Multiplier      Final Value
  -------   ---------------------------      -----------
  H-H       (15 + 15) * 1.0                  30
  C-H       (35 + 15) * 1.0                  50  (baseline)
  N-H       (35 + 15) * 1.0                  50
  O-H       (35 + 15) * 1.0                  50
  P-H       (40 + 15) * 1.0                  55

  C-C       (35 + 35) * 1.0                  70  (baseline)
  C=C       (35 + 35) * 0.9                  63
  C≡C       (35 + 35) * 0.85                 60
  
  C-O       (35 + 35) * 1.0                  70
  C=O       (35 + 35) * 0.9                  63

  S-O       (40 + 35) * 1.0                  75
  S=O       (40 + 35) * 0.9                  68  (SO4, SO3)
  
  P-Cl      (40 + 40) * 1.0                  80  (PCl3)
  Xe=O      (50 + 35) * 0.9                  76  (XeO3)
  
  F-F       (35 + 35) * 1.0                  70
  Cl-Cl     (40 + 40) * 1.0                  80
  I-I       (50 + 50) * 1.0                  100
 ==========================================================================
*/


const ELEMENT_PROPS = {
    // --- Period 1 ---
    "H":  { ve: 1, c3d: "#F0F0F0", r3d: 12, lp: 0, mass: 1.008, en: 2.20 }, // [保留] 白灰
    "He": { ve: 2, c3d: "#A5F3FC", r3d: 11, lp: 1, mass: 4.002, en: 0 },    // [保留] 淡青

    // --- Period 2 (radius decreases: Li > Be > B > C > N > O > F) ---
    "Li": { ve: 1, c3d: "#E879F9", r3d: 26, lp: 0, mass: 6.94, en: 0.98 },  // [保留] 粉紫
    "Be": { ve: 2, c3d: "#C2F970", r3d: 22, lp: 0, mass: 9.012, en: 1.57 }, // [保留] 萊姆綠
    "B":  { ve: 3, c3d: "#FDBA74", r3d: 20, lp: 0, mass: 10.81, en: 2.04 }, // [保留] 蜜桃橘 (您喜歡的顏色)
    "C":  { ve: 4, c3d: "#94A3B8", r3d: 19, lp: 0, mass: 12.011, en: 2.55 },// [reserved] blue-gray
    "N":  { ve: 5, c3d: "#3B82F6", r3d: 18, lp: 1, mass: 14.007, en: 3.04 },// [reserved] bright-blue
    "O":  { ve: 6, c3d: "#EF4444", r3d: 17, lp: 2, mass: 15.999, en: 3.44 },// [reserved] red
    "F":  { ve: 7, c3d: "#90E050", r3d: 16, lp: 3, mass: 18.998, en: 3.98 },// [reserved] fresh-green
    "Ne": { ve: 8, c3d: "#67E8F9", r3d: 15, lp: 4, mass: 20.180, en: 0 },   // [保留] 青

    // --- Period 3 (radius decreases: Na > Mg > Al > Si > P > S > Cl) ---
    "Na": { ve: 1, c3d: "#C084FC", r3d: 30, lp: 0, mass: 22.990, en: 0.93 },// [保留] 紫
    "Mg": { ve: 2, c3d: "#10B981", r3d: 26, lp: 0, mass: 24.305, en: 1.31 },// [保留] 翡翠綠
    "Al": { ve: 3, c3d: "#E2E8F0", r3d: 24, lp: 0, mass: 26.982, en: 1.61 },// [保留] 淺灰
    "Si": { ve: 4, c3d: "#CBD5E1", r3d: 22, lp: 0, mass: 28.085, en: 1.90 },// [保留] 灰
    "P":  { ve: 5, c3d: "#F97316", r3d: 21, lp: 0, mass: 30.974, en: 2.19 },// [保留] 橘
    "S":  { ve: 6, c3d: "#FACC15", r3d: 20, lp: 0, mass: 32.06, en: 2.58 }, // [保留] 黃
    "Cl": { ve: 7, c3d: "#22C55E", r3d: 19, lp: 3, mass: 35.45, en: 3.16 }, // [保留] 深綠 (比F深，符合視覺邏輯)
    "Ar": { ve: 8, c3d: "#38BDF8", r3d: 18, lp: 4, mass: 39.948, en: 0 },   // [保留] 天藍

    // --- Period 4 (K > Ca > Sc ... > Br) ---
    "K":  { ve: 1, c3d: "#8B5CF6", r3d: 36, lp: 0, mass: 39.098, en: 0.82 },// [保留] 靛紫
    "Ca": { ve: 2, c3d: "#4ADE80", r3d: 32, lp: 0, mass: 40.078, en: 1.00 },// [保留] 淺綠
    "Sc": { ve: 3, c3d: "#E6E6E6", r3d: 28, lp: 0, mass: 44.96, en: 1.36 }, // [added] CPK silver-white
    "Ti": { ve: 4, c3d: "#BFC2C7", r3d: 26, lp: 0, mass: 47.87, en: 1.54 }, // [added] CPK titanium-gray
    "V":  { ve: 5, c3d: "#A6A6AB", r3d: 25, lp: 0, mass: 50.94, en: 1.63 }, // [新增] 灰
    "Cr": { ve: 6, c3d: "#94A3B8", r3d: 24, lp: 0, mass: 51.996, en: 1.66 },// [保留] 鉻灰
    "Mn": { ve: 7, c3d: "#D946EF", r3d: 24, lp: 0, mass: 54.938, en: 1.55 },// [保留] 紫紅 (很適合Mn)
    "Fe": { ve: 8, c3d: "#EA580C", r3d: 24, lp: 0, mass: 55.845, en: 1.83 },// [保留] 鐵鏽橘
    "Co": { ve: 9, c3d: "#F472B6", r3d: 23, lp: 0, mass: 58.93, en: 1.88 }, // [新增] CPK 粉紅 (鈷)
    "Ni": { ve: 10, c3d: "#50D050", r3d: 23, lp: 0, mass: 58.69, en: 1.91 },// [新增] CPK 綠 (鎳)
    "Cu": { ve: 11, c3d: "#D97706", r3d: 23, lp: 0, mass: 63.546, en: 1.90 },// [保留] 銅色
    "Zn": { ve: 12, c3d: "#78716C", r3d: 23, lp: 0, mass: 65.38, en: 1.65 },// [保留] 鋅灰
    "Ga": { ve: 3, c3d: "#C28F8F", r3d: 23, lp: 0, mass: 69.72, en: 1.81 }, // [新增] 紅褐
    "Ge": { ve: 4, c3d: "#668F8F", r3d: 23, lp: 0, mass: 72.63, en: 2.01 }, // [新增] 灰綠
    "As": { ve: 5, c3d: "#BD80E3", r3d: 22, lp: 0, mass: 74.92, en: 2.18 }, // [updated] purple (distinct from Na)
    "Se": { ve: 6, c3d: "#FFA100", r3d: 22, lp: 0, mass: 78.96, en: 2.55 }, // [updated] deep-orange (distinct from P)
    "Br": { ve: 7, c3d: "#B91C1C", r3d: 21, lp: 3, mass: 79.904, en: 2.96 },// [保留] 深紅
    "Kr": { ve: 8, c3d: "#5CB8D1", r3d: 20, lp: 4, mass: 83.80, en: 3.00 }, // [新增] 青

    // --- Other common elements (Periods 5, 6) ---
    "Rb": { ve: 1, c3d: "#702EB0", r3d: 38, lp: 0, mass: 85.47, en: 0.82 }, // CPK 紫
    "Sr": { ve: 2, c3d: "#00FF00", r3d: 34, lp: 0, mass: 87.62, en: 0.95 }, // CPK 綠
    "Ag": { ve: 11, c3d: "#F1F5F9", r3d: 25, lp: 0, mass: 107.87, en: 1.93 },// [保留] 銀白
    "Sn": { ve: 4, c3d: "#668080", r3d: 25, lp: 0, mass: 118.7, en: 1.96 }, // 灰
    "Sb": { ve: 5, c3d: "#A855F7", r3d: 25, lp: 0, mass: 121.76, en: 2.05 },// [保留] 紫
    "Te": { ve: 6, c3d: "#EA580C", r3d: 25, lp: 0, mass: 127.60, en: 2.10 },// [保留] 橘褐
    "I":  { ve: 7, c3d: "#A855F7", r3d: 24, lp: 3, mass: 126.90, en: 2.66 },// [保留] 紫
    "Xe": { ve: 8, c3d: "#818CF8", r3d: 24, lp: 3, mass: 131.29, en: 2.60 },// [保留] 藍紫
    "Cs": { ve: 1, c3d: "#57178F", r3d: 42, lp: 0, mass: 132.9, en: 0.79 }, // CPK 深紫
    "Ba": { ve: 2, c3d: "#00C900", r3d: 38, lp: 0, mass: 137.3, en: 0.89 }, // CPK 深綠
    "Pt": { ve: 10, c3d: "#D0D0E0", r3d: 25, lp: 0, mass: 195.1, en: 2.28 },// 鉑
    "Au": { ve: 11, c3d: "#F59E0B", r3d: 25, lp: 0, mass: 196.97, en: 2.54 },// [保留] 金黃
    "Hg": { ve: 12, c3d: "#B8B8D0", r3d: 24, lp: 0, mass: 200.6, en: 2.00 },// 汞
    "Pb": { ve: 4, c3d: "#575961", r3d: 26, lp: 0, mass: 207.2, en: 2.33 }, // 鉛

    // --- [Supplement] Metals and radioactive elements (for crystal structures) ---
    "Po": { ve: 6, c3d: "#AB5C00", r3d: 26, lp: 0, mass: 209, en: 2.0 },    // [新增] 釙 (金屬) - 深橘褐
    "Fr": { ve: 1, c3d: "#420066", r3d: 44, lp: 0, mass: 223, en: 0.7 },    // [新增] 鍅 (1A) - 極深紫
    "Ra": { ve: 2, c3d: "#006400", r3d: 40, lp: 0, mass: 226, en: 0.9 },    // [新增] 鐳 (2A) - 深綠
// --- Materials for electron orbitals (no text) ---
    // Tip: use different amounts of spaces as IDs so no text appears, but colors differ
    " ":      { ve: 0, c3d: "#3B82F6", r3d: 0 }, // s (藍)
    "  ":     { ve: 0, c3d: "#10B981", r3d: 0 }, // p (備用)
    "   ":    { ve: 0, c3d: "#F59E0B", r3d: 0 }, // d (橘)
    
    // --- Coordinate axes system ---
    "Origin": { ve: 0, c3d: "#000000", r3d: 0,    lp: 0, mass: 0, en: 0 }, // 隱藏原點
    "Axis":   { ve: 0, c3d: "#444444", r3d: 1.0,  lp: 0, mass: 0, en: 0 }, // 極細深灰軸
    "X":      { ve: 0, c3d: "#EF4444", r3d: 0,    lp: 0, mass: 0, en: 0 }, // 紅色 X
    "Y":      { ve: 0, c3d: "#22C55E", r3d: 0,    lp: 0, mass: 0, en: 0 }, // 綠色 Y
    "Z":      { ve: 0, c3d: "#3B82F6", r3d: 0,    lp: 0, mass: 0, en: 0 }  // 藍色 Z
};
// hide origin, axes remain for reference

// ========== [Insert this JS] Electron configuration data ==========
    const ELECTRON_DATA = [
    // Period 1
    { z: 1, s: "H", n: "Hydrogen", cn: "Hydrogen", type: "nonmetal", state: "gas", mp: "-259°C", bp: "-253°C", p: 1, g: "1A", iupac: 1, c: "1s1", noble: "1s1" },
    { z: 2, s: "He", n: "Helium", cn: "Helium", type: "nonmetal", state: "gas", mp: "-272°C", bp: "-269°C", p: 1, g: "8A", iupac: 18, c: "1s2", noble: "1s2" },
    // Period 2
    { z: 3, s: "Li", n: "Lithium", cn: "Lithium", type: "metal", state: "solid", mp: "180°C", bp: "1342°C", p: 2, g: "1A", iupac: 1, c: "1s2 2s1", noble: "[He] 2s1" },
    { z: 4, s: "Be", n: "Beryllium", cn: "Beryllium", type: "metal", state: "solid", mp: "1287°C", bp: "2469°C", p: 2, g: "2A", iupac: 2, c: "1s2 2s2", noble: "[He] 2s2" },
    { z: 5, s: "B", n: "Boron", cn: "Boron", type: "metalloid", state: "solid", mp: "2076°C", bp: "3927°C", p: 2, g: "3A", iupac: 13, c: "1s2 2s2 2p1", noble: "[He] 2s2 2p1" },
    { z: 6, s: "C", n: "Carbon", cn: "Carbon", type: "nonmetal", state: "solid", mp: "3550°C", bp: "4027°C", p: 2, g: "4A", iupac: 14, c: "1s2 2s2 2p2", noble: "[He] 2s2 2p2" },
    { z: 7, s: "N", n: "Nitrogen", cn: "Nitrogen", type: "nonmetal", state: "gas", mp: "-210°C", bp: "-196°C", p: 2, g: "5A", iupac: 15, c: "1s2 2s2 2p3", noble: "[He] 2s2 2p3" },
    { z: 8, s: "O", n: "Oxygen", cn: "Oxygen", type: "nonmetal", state: "gas", mp: "-218°C", bp: "-183°C", p: 2, g: "6A", iupac: 16, c: "1s2 2s2 2p4", noble: "[He] 2s2 2p4" },
    { z: 9, s: "F", n: "Fluorine", cn: "Fluorine", type: "nonmetal", state: "gas", mp: "-220°C", bp: "-188°C", p: 2, g: "7A", iupac: 17, c: "1s2 2s2 2p5", noble: "[He] 2s2 2p5" },
    { z: 10, s: "Ne", n: "Neon", cn: "Neon", type: "nonmetal", state: "gas", mp: "-249°C", bp: "-246°C", p: 2, g: "8A", iupac: 18, c: "1s2 2s2 2p6", noble: "[He] 2s2 2p6" },
    // Period 3
    { z: 11, s: "Na", n: "Sodium", cn: "Sodium", type: "metal", state: "solid", mp: "98°C", bp: "883°C", p: 3, g: "1A", iupac: 1, c: "1s2 2s2 2p6 3s1", noble: "[Ne] 3s1" },
    { z: 12, s: "Mg", n: "Magnesium", cn: "Magnesium", type: "metal", state: "solid", mp: "650°C", bp: "1090°C", p: 3, g: "2A", iupac: 2, c: "1s2 2s2 2p6 3s2", noble: "[Ne] 3s2" },
    { z: 13, s: "Al", n: "Aluminium", cn: "Aluminium", type: "metal", state: "solid", mp: "660°C", bp: "2519°C", p: 3, g: "3A", iupac: 13, c: "1s2 2s2 2p6 3s2 3p1", noble: "[Ne] 3s2 3p1" },
    { z: 14, s: "Si", n: "Silicon", cn: "Silicon", type: "metalloid", state: "solid", mp: "1414°C", bp: "3265°C", p: 3, g: "4A", iupac: 14, c: "1s2 2s2 2p6 3s2 3p2", noble: "[Ne] 3s2 3p2" },
    { z: 15, s: "P", n: "Phosphorus", cn: "Phosphorus", type: "nonmetal", state: "solid", mp: "44°C", bp: "280°C", p: 3, g: "5A", iupac: 15, c: "1s2 2s2 2p6 3s2 3p3", noble: "[Ne] 3s2 3p3" },
    { z: 16, s: "S", n: "Sulfur", cn: "Sulfur", type: "nonmetal", state: "solid", mp: "115°C", bp: "445°C", p: 3, g: "6A", iupac: 16, c: "1s2 2s2 2p6 3s2 3p4", noble: "[Ne] 3s2 3p4" },
    { z: 17, s: "Cl", n: "Chlorine", cn: "Chlorine", type: "nonmetal", state: "gas", mp: "-101°C", bp: "-34°C", p: 3, g: "7A", iupac: 17, c: "1s2 2s2 2p6 3s2 3p5", noble: "[Ne] 3s2 3p5" },
    { z: 18, s: "Ar", n: "Argon", cn: "Argon", type: "nonmetal", state: "gas", mp: "-189°C", bp: "-186°C", p: 3, g: "8A", iupac: 18, c: "1s2 2s2 2p6 3s2 3p6", noble: "[Ne] 3s2 3p6" },
    // Period 4
    { z: 19, s: "K", n: "Potassium", cn: "Potassium", type: "metal", state: "solid", mp: "63°C", bp: "759°C", p: 4, g: "1A", iupac: 1, c: "1s2 2s2 2p6 3s2 3p6 4s1", noble: "[Ar] 4s1" },
    { z: 20, s: "Ca", n: "Calcium", cn: "Calcium", type: "metal", state: "solid", mp: "842°C", bp: "1484°C", p: 4, g: "2A", iupac: 2, c: "1s2 2s2 2p6 3s2 3p6 4s2", noble: "[Ar] 4s2" },
    { z: 21, s: "Sc", n: "Scandium", cn: "Scandium", type: "metal", state: "solid", mp: "1541°C", bp: "2836°C", p: 4, g: "3B", iupac: 3, c: "1s2 2s2 2p6 3s2 3p6 3d1 4s2", noble: "[Ar] 3d1 4s2" },
    { z: 22, s: "Ti", n: "Titanium", cn: "Titanium", type: "metal", state: "solid", mp: "1668°C", bp: "3287°C", p: 4, g: "4B", iupac: 4, c: "1s2 2s2 2p6 3s2 3p6 3d2 4s2", noble: "[Ar] 3d2 4s2" },
    { z: 23, s: "V", n: "Vanadium", cn: "Vanadium", type: "metal", state: "solid", mp: "1910°C", bp: "3407°C", p: 4, g: "5B", iupac: 5, c: "1s2 2s2 2p6 3s2 3p6 3d3 4s2", noble: "[Ar] 3d3 4s2" },
    { z: 24, s: "Cr", n: "Chromium", cn: "Chromium", type: "metal", state: "solid", mp: "1907°C", bp: "2671°C", p: 4, g: "6B", iupac: 6, c: "1s2 2s2 2p6 3s2 3p6 3d5 4s1", noble: "[Ar] 3d5 4s1", ex: true },    
    { z: 25, s: "Mn", n: "Manganese", cn: "Manganese", type: "metal", state: "solid", mp: "1246°C", bp: "2061°C", p: 4, g: "7B", iupac: 7, c: "1s2 2s2 2p6 3s2 3p6 3d5 4s2", noble: "[Ar] 3d5 4s2" },
    { z: 26, s: "Fe", n: "Iron", cn: "Iron", type: "metal", state: "solid", mp: "1538°C", bp: "2861°C", p: 4, g: "8B", iupac: 8, c: "1s2 2s2 2p6 3s2 3p6 3d6 4s2", noble: "[Ar] 3d6 4s2" },
    { z: 27, s: "Co", n: "Cobalt", cn: "Cobalt", type: "metal", state: "solid", mp: "1495°C", bp: "2927°C", p: 4, g: "8B", iupac: 9, c: "1s2 2s2 2p6 3s2 3p6 3d7 4s2", noble: "[Ar] 3d7 4s2" },
    { z: 28, s: "Ni", n: "Nickel", cn: "Nickel", type: "metal", state: "solid", mp: "1455°C", bp: "2730°C", p: 4, g: "8B", iupac: 10, c: "1s2 2s2 2p6 3s2 3p6 3d8 4s2", noble: "[Ar] 3d8 4s2" },
    { z: 29, s: "Cu", n: "Copper", cn: "Copper", type: "metal", state: "solid", mp: "1085°C", bp: "2562°C", p: 4, g: "1B", iupac: 11, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s1", noble: "[Ar] 3d10 4s1", ex: true },
    { z: 30, s: "Zn", n: "Zinc", cn: "Zinc", type: "metal", state: "solid", mp: "420°C", bp: "907°C", p: 4, g: "2B", iupac: 12, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2", noble: "[Ar] 3d10 4s2" },
    { z: 31, s: "Ga", n: "Gallium", cn: "Gallium", type: "metal", state: "solid", mp: "30°C", bp: "2204°C", p: 4, g: "3A", iupac: 13, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p1", noble: "[Ar] 3d10 4s2 4p1" },
    { z: 32, s: "Ge", n: "Germanium", cn: "Germanium", type: "metalloid", state: "solid", mp: "938°C", bp: "2833°C", p: 4, g: "4A", iupac: 14, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p2", noble: "[Ar] 3d10 4s2 4p2" },
    { z: 33, s: "As", n: "Arsenic", cn: "Arsenic", type: "metalloid", state: "solid", mp: "817°C", bp: "614°C", p: 4, g: "5A", iupac: 15, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p3", noble: "[Ar] 3d10 4s2 4p3" },
    { z: 34, s: "Se", n: "Selenium", cn: "Selenium", type: "nonmetal", state: "solid", mp: "221°C", bp: "685°C", p: 4, g: "6A", iupac: 16, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p4", noble: "[Ar] 3d10 4s2 4p4" },
    { z: 35, s: "Br", n: "Bromine", cn: "Bromine", type: "nonmetal", state: "liquid", mp: "-7°C", bp: "59°C", p: 4, g: "7A", iupac: 17, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p5", noble: "[Ar] 3d10 4s2 4p5" },
    { z: 36, s: "Kr", n: "Krypton", cn: "Krypton", type: "nonmetal", state: "gas", mp: "-157°C", bp: "-153°C", p: 4, g: "8A", iupac: 18, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6", noble: "[Ar] 3d10 4s2 4p6" },
    // Period 5
    { z: 37, s: "Rb", n: "Rubidium", cn: "Rubidium", type: "metal", state: "solid", mp: "39°C", bp: "688°C", p: 5, g: "1A", iupac: 1, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 5s1", noble: "[Kr] 5s1" },
    { z: 38, s: "Sr", n: "Strontium", cn: "Strontium", type: "metal", state: "solid", mp: "777°C", bp: "1382°C", p: 5, g: "2A", iupac: 2, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 5s2", noble: "[Kr] 5s2" },
    { z: 39, s: "Y", n: "Yttrium", cn: "Yttrium", type: "metal", state: "solid", mp: "1526°C", bp: "3338°C", p: 5, g: "3B", iupac: 3, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d1 5s2", noble: "[Kr] 4d1 5s2" },
    { z: 40, s: "Zr", n: "Zirconium", cn: "Zirconium", type: "metal", state: "solid", mp: "1855°C", bp: "4409°C", p: 5, g: "4B", iupac: 4, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d2 5s2", noble: "[Kr] 4d2 5s2" },
    { z: 41, s: "Nb", n: "Niobium", cn: "Niobium", type: "metal", state: "solid", mp: "2477°C", bp: "4744°C", p: 5, g: "5B", iupac: 5, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d4 5s1", noble: "[Kr] 4d4 5s1", ex: true },    
    { z: 42, s: "Mo", n: "Molybdenum", cn: "Molybdenum", type: "metal", state: "solid", mp: "2623°C", bp: "4639°C", p: 5, g: "6B", iupac: 6, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d5 5s1", noble: "[Kr] 4d5 5s1", ex: true },    
    { z: 43, s: "Tc", n: "Technetium", cn: "Technetium", type: "metal", state: "solid", mp: "2157°C", bp: "4265°C", p: 5, g: "7B", iupac: 7, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d5 5s2", noble: "[Kr] 4d5 5s2" },
    { z: 44, s: "Ru", n: "Ruthenium", cn: "Ruthenium", type: "metal", state: "solid", mp: "2334°C", bp: "4150°C", p: 5, g: "8B", iupac: 8, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d7 5s1", noble: "[Kr] 4d7 5s1" },
    { z: 45, s: "Rh", n: "Rhodium", cn: "Rhodium", type: "metal", state: "solid", mp: "1964°C", bp: "3695°C", p: 5, g: "8B", iupac: 9, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d8 5s1", noble: "[Kr] 4d8 5s1", ex: true },    
    { z: 46, s: "Pd", n: "Palladium", cn: "Palladium", type: "metal", state: "solid", mp: "1555°C", bp: "2963°C", p: 5, g: "8B", iupac: 10, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10", noble: "[Kr] 4d10", ex: true },    
    { z: 47, s: "Ag", n: "Silver", cn: "Silver", type: "metal", state: "solid", mp: "962°C", bp: "2162°C", p: 5, g: "1B", iupac: 11, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s1", noble: "[Kr] 4d10 5s1", ex: true },    
    { z: 48, s: "Cd", n: "Cadmium", cn: "Cadmium", type: "metal", state: "solid", mp: "321°C", bp: "767°C", p: 5, g: "2B", iupac: 12, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2", noble: "[Kr] 4d10 5s2" },
    { z: 49, s: "In", n: "Indium", cn: "Indium", type: "metal", state: "solid", mp: "157°C", bp: "2072°C", p: 5, g: "3A", iupac: 13, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p1", noble: "[Kr] 4d10 5s2 5p1" },
    { z: 50, s: "Sn", n: "Tin", cn: "Tin", type: "metal", state: "solid", mp: "232°C", bp: "2602°C", p: 5, g: "4A", iupac: 14, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p2", noble: "[Kr] 4d10 5s2 5p2" },
    { z: 51, s: "Sb", n: "Antimony", cn: "Antimony", type: "metalloid", state: "solid", mp: "631°C", bp: "1587°C", p: 5, g: "5A", iupac: 15, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p3", noble: "[Kr] 4d10 5s2 5p3" },
    { z: 52, s: "Te", n: "Tellurium", cn: "Tellurium", type: "metalloid", state: "solid", mp: "450°C", bp: "988°C", p: 5, g: "6A", iupac: 16, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p4", noble: "[Kr] 4d10 5s2 5p4" },
    { z: 53, s: "I", n: "Iodine", cn: "Iodine", type: "nonmetal", state: "solid", mp: "114°C", bp: "184°C", p: 5, g: "7A", iupac: 17, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p5", noble: "[Kr] 4d10 5s2 5p5" },
    { z: 54, s: "Xe", n: "Xenon", cn: "Xenon", type: "nonmetal", state: "gas", mp: "-112°C", bp: "-108°C", p: 5, g: "8A", iupac: 18, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6", noble: "[Kr] 4d10 5s2 5p6" },
    // Period 6
    { z: 55, s: "Cs", n: "Cesium", cn: "Cesium", type: "metal", state: "solid", mp: "28°C", bp: "671°C", p: 6, g: "1A", iupac: 1, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 6s1", noble: "[Xe] 6s1" },
    { z: 56, s: "Ba", n: "Barium", cn: "Barium", type: "metal", state: "solid", mp: "727°C", bp: "1897°C", p: 6, g: "2A", iupac: 2, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 6s2", noble: "[Xe] 6s2" },
    { z: 57, s: "La", n: "Lanthanum", cn: "Lanthanum", type: "metal", state: "solid", mp: "920°C", bp: "3464°C", p: 6, g: "lanthanide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 5d1 6s2", noble: "[Xe] 5d1 6s2", ex: true },    
    { z: 58, s: "Ce", n: "Cerium", cn: "Cerium", type: "metal", state: "solid", mp: "795°C", bp: "3443°C", p: 6, g: "lanthanide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f1 5d1 6s2", noble: "[Xe] 4f1 5d1 6s2", ex: true },    
    { z: 59, s: "Pr", n: "Praseodymium", cn: "Praseodymium", type: "metal", state: "solid", mp: "931°C", bp: "3520°C", p: 6, g: "lanthanide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f3 6s2", noble: "[Xe] 4f3 6s2" },
    { z: 60, s: "Nd", n: "Neodymium", cn: "Neodymium", type: "metal", state: "solid", mp: "1024°C", bp: "3074°C", p: 6, g: "lanthanide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f4 6s2", noble: "[Xe] 4f4 6s2" },
    { z: 61, s: "Pm", n: "Promethium", cn: "Promethium", type: "metal", state: "solid", mp: "1042°C", bp: "3000°C", p: 6, g: "lanthanide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f5 6s2", noble: "[Xe] 4f5 6s2" },
    { z: 62, s: "Sm", n: "Samarium", cn: "Samarium", type: "metal", state: "solid", mp: "1072°C", bp: "1794°C", p: 6, g: "lanthanide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f6 6s2", noble: "[Xe] 4f6 6s2" },
    { z: 63, s: "Eu", n: "Europium", cn: "Europium", type: "metal", state: "solid", mp: "826°C", bp: "1529°C", p: 6, g: "lanthanide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f7 6s2", noble: "[Xe] 4f7 6s2" },
    { z: 64, s: "Gd", n: "Gadolinium", cn: "Gadolinium", type: "metal", state: "solid", mp: "1312°C", bp: "3273°C", p: 6, g: "lanthanide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f7 5d1 6s2", noble: "[Xe] 4f7 5d1 6s2", ex: true },    
    { z: 65, s: "Tb", n: "Terbium", cn: "Terbium", type: "metal", state: "solid", mp: "1356°C", bp: "3230°C", p: 6, g: "lanthanide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f9 6s2", noble: "[Xe] 4f9 6s2" },
    { z: 66, s: "Dy", n: "Dysprosium", cn: "Dysprosium", type: "metal", state: "solid", mp: "1407°C", bp: "2567°C", p: 6, g: "lanthanide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f10 6s2", noble: "[Xe] 4f10 6s2" },
    { z: 67, s: "Ho", n: "Holmium", cn: "Holmium", type: "metal", state: "solid", mp: "1461°C", bp: "2720°C", p: 6, g: "lanthanide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f11 6s2", noble: "[Xe] 4f11 6s2" },
    { z: 68, s: "Er", n: "Erbium", cn: "Erbium", type: "metal", state: "solid", mp: "1529°C", bp: "2868°C", p: 6, g: "lanthanide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f12 6s2", noble: "[Xe] 4f12 6s2" },
    { z: 69, s: "Tm", n: "Thulium", cn: "Thulium", type: "metal", state: "solid", mp: "1545°C", bp: "1950°C", p: 6, g: "lanthanide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f13 6s2", noble: "[Xe] 4f13 6s2" },
    { z: 70, s: "Yb", n: "Ytterbium", cn: "Ytterbium", type: "metal", state: "solid", mp: "824°C", bp: "1196°C", p: 6, g: "lanthanide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 6s2", noble: "[Xe] 4f14 6s2" },
    { z: 71, s: "Lu", n: "Lutetium", cn: "Lutetium", type: "metal", state: "solid", mp: "1663°C", bp: "3402°C", p: 6, g: "lanthanide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d1 6s2", noble: "[Xe] 4f14 5d1 6s2" },
    { z: 72, s: "Hf", n: "Hafnium", cn: "Hafnium", type: "metal", state: "solid", mp: "2233°C", bp: "4603°C", p: 6, g: "4B", iupac: 4, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d2 6s2", noble: "[Xe] 4f14 5d2 6s2" },
    { z: 73, s: "Ta", n: "Tantalum", cn: "Tantalum", type: "metal", state: "solid", mp: "3017°C", bp: "5458°C", p: 6, g: "5B", iupac: 5, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d3 6s2", noble: "[Xe] 4f14 5d3 6s2" },
    { z: 74, s: "W", n: "Tungsten", cn: "Tungsten", type: "metal", state: "solid", mp: "3422°C", bp: "5930°C", p: 6, g: "6B", iupac: 6, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d4 6s2", noble: "[Xe] 4f14 5d4 6s2" },
    { z: 75, s: "Re", n: "Rhenium", cn: "Rhenium", type: "metal", state: "solid", mp: "3186°C", bp: "5596°C", p: 6, g: "7B", iupac: 7, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d5 6s2", noble: "[Xe] 4f14 5d5 6s2" },
    { z: 76, s: "Os", n: "Osmium", cn: "Osmium", type: "metal", state: "solid", mp: "3033°C", bp: "5012°C", p: 6, g: "8B", iupac: 8, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d6 6s2", noble: "[Xe] 4f14 5d6 6s2" },
    { z: 77, s: "Ir", n: "Iridium", cn: "Iridium", type: "metal", state: "solid", mp: "2446°C", bp: "4428°C", p: 6, g: "8B", iupac: 9, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d7 6s2", noble: "[Xe] 4f14 5d7 6s2" },
    { z: 78, s: "Pt", n: "Platinum", cn: "Platinum", type: "metal", state: "solid", mp: "1768°C", bp: "3825°C", p: 6, g: "8B", iupac: 10, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d9 6s1", noble: "[Xe] 4f14 5d9 6s1", ex: true },    
    { z: 79, s: "Au", n: "Gold", cn: "Gold", type: "metal", state: "solid", mp: "1064°C", bp: "2970°C", p: 6, g: "1B", iupac: 11, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s1", noble: "[Xe] 4f14 5d10 6s1", ex: true },    
    { z: 80, s: "Hg", n: "Mercury", cn: "Mercury", type: "metal", state: "liquid", mp: "-39°C", bp: "357°C", p: 6, g: "2B", iupac: 12, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2", noble: "[Xe] 4f14 5d10 6s2" },
    { z: 81, s: "Tl", n: "Thallium", cn: "Thallium", type: "metal", state: "solid", mp: "304°C", bp: "1473°C", p: 6, g: "3A", iupac: 13, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p1", noble: "[Xe] 4f14 5d10 6s2 6p1" },
    { z: 82, s: "Pb", n: "Lead", cn: "Lead", type: "metal", state: "solid", mp: "327°C", bp: "1749°C", p: 6, g: "4A", iupac: 14, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p2", noble: "[Xe] 4f14 5d10 6s2 6p2" },
    { z: 83, s: "Bi", n: "Bismuth", cn: "Bismuth", type: "metal", state: "solid", mp: "271°C", bp: "1564°C", p: 6, g: "5A", iupac: 15, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p3", noble: "[Xe] 4f14 5d10 6s2 6p3" },
    { z: 84, s: "Po", n: "Polonium", cn: "Polonium", type: "metal", state: "solid", mp: "254°C", bp: "962°C", p: 6, g: "6A", iupac: 16, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p4", noble: "[Xe] 4f14 5d10 6s2 6p4" },
    { z: 85, s: "At", n: "Astatine", cn: "Astatine", type: "metalloid", state: "solid", mp: "302°C", bp: "337°C", p: 6, g: "7A", iupac: 17, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p5", noble: "[Xe] 4f14 5d10 6s2 6p5" },
    { z: 86, s: "Rn", n: "Radon", cn: "Radon", type: "nonmetal", state: "gas", mp: "-71°C", bp: "-62°C", p: 6, g: "8A", iupac: 18, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6", noble: "[Xe] 4f14 5d10 6s2 6p6" },
    // Period 7
    { z: 87, s: "Fr", n: "Francium", cn: "Francium", type: "metal", state: "solid", mp: "27°C", bp: "677°C", p: 7, g: "1A", iupac: 1, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 7s1", noble: "[Rn] 7s1" },
    { z: 88, s: "Ra", n: "Radium", cn: "Radium", type: "metal", state: "solid", mp: "700°C", bp: "1737°C", p: 7, g: "2A", iupac: 2, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 7s2", noble: "[Rn] 7s2" },
    { z: 89, s: "Ac", n: "Actinium", cn: "Actinium", type: "metal", state: "solid", mp: "1050°C", bp: "3198°C", p: 7, g: "actinide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 6d1 7s2", noble: "[Rn] 6d1 7s2", ex: true },
    { z: 90, s: "Th", n: "Thorium", cn: "Thorium", type: "metal", state: "solid", mp: "1750°C", bp: "4788°C", p: 7, g: "actinide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 6d2 7s2", noble: "[Rn] 6d2 7s2", ex: true },
    { z: 91, s: "Pa", n: "Protactinium", cn: "Protactinium", type: "metal", state: "solid", mp: "1568°C", bp: "4027°C", p: 7, g: "actinide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f2 6d1 7s2", noble: "[Rn] 5f2 6d1 7s2", ex: true },
    { z: 92, s: "U", n: "Uranium", cn: "Uranium", type: "metal", state: "solid", mp: "1132°C", bp: "4131°C", p: 7, g: "actinide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f3 6d1 7s2", noble: "[Rn] 5f3 6d1 7s2", ex: true },
    { z: 93, s: "Np", n: "Neptunium", cn: "Neptunium", type: "metal", state: "solid", mp: "644°C", bp: "3902°C", p: 7, g: "actinide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f4 6d1 7s2", noble: "[Rn] 5f4 6d1 7s2", ex: true },
    { z: 94, s: "Pu", n: "Plutonium", cn: "Plutonium", type: "metal", state: "solid", mp: "640°C", bp: "3228°C", p: 7, g: "actinide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f6 7s2", noble: "[Rn] 5f6 7s2" },
    { z: 95, s: "Am", n: "Americium", cn: "Americium", type: "metal", state: "solid", mp: "1176°C", bp: "2607°C", p: 7, g: "actinide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f7 7s2", noble: "[Rn] 5f7 7s2" },
    { z: 96, s: "Cm", n: "Curium", cn: "Curium", type: "metal", state: "solid", mp: "1340°C", bp: "3110°C", p: 7, g: "actinide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f7 6d1 7s2", noble: "[Rn] 5f7 6d1 7s2", ex: true },    
    { z: 97, s: "Bk", n: "Berkelium", cn: "Berkelium", type: "metal", state: "solid", mp: "986°C", bp: "-", p: 7, g: "actinide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f9 7s2", noble: "[Rn] 5f9 7s2" },
    { z: 98, s: "Cf", n: "Californium", cn: "Californium", type: "metal", state: "solid", mp: "900°C", bp: "-", p: 7, g: "actinide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f10 7s2", noble: "[Rn] 5f10 7s2" },
    { z: 99, s: "Es", n: "Einsteinium", cn: "Einsteinium", type: "metal", state: "solid", mp: "860°C", bp: "-", p: 7, g: "actinide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f11 7s2", noble: "[Rn] 5f11 7s2" },
    { z: 100, s: "Fm", n: "Fermium", cn: "Fermium", type: "metal", state: "solid", mp: "1527°C", bp: "-", p: 7, g: "actinide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f12 7s2", noble: "[Rn] 5f12 7s2" },
    { z: 101, s: "Md", n: "Mendelevium", cn: "Mendelevium", type: "metal", state: "solid", mp: "827°C", bp: "-", p: 7, g: "actinide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f13 7s2", noble: "[Rn] 5f13 7s2" },
    { z: 102, s: "No", n: "Nobelium", cn: "Nobelium", type: "metal", state: "solid", mp: "827°C", bp: "-", p: 7, g: "actinide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f14 7s2", noble: "[Rn] 5f14 7s2" },
    { z: 103, s: "Lr", n: "Lawrencium", cn: "Lawrencium", type: "metal", state: "solid", mp: "1627°C", bp: "-", p: 7, g: "actinide", iupac: "-", c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f14 7s2 7p1", noble: "[Rn] 5f14 7s2 7p1" },
    { z: 104, s: "Rf", n: "Rutherfordium", cn: "Rutherfordium", type: "metal", state: "solid", mp: "-", bp: "-", p: 7, g: "4B", iupac: 4, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f14 6d2 7s2", noble: "[Rn] 5f14 6d2 7s2" },
    { z: 105, s: "Db", n: "Dubnium", cn: "Dubnium", type: "metal", state: "solid", mp: "-", bp: "-", p: 7, g: "5B", iupac: 5, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f14 6d3 7s2", noble: "[Rn] 5f14 6d3 7s2" },
    { z: 106, s: "Sg", n: "Seaborgium", cn: "Seaborgium", type: "metal", state: "solid", mp: "-", bp: "-", p: 7, g: "6B", iupac: 6, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f14 6d4 7s2", noble: "[Rn] 5f14 6d4 7s2" },
    { z: 107, s: "Bh", n: "Bohrium", cn: "Bohrium", type: "metal", state: "solid", mp: "-", bp: "-", p: 7, g: "7B", iupac: 7, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f14 6d5 7s2", noble: "[Rn] 5f14 6d5 7s2" },
    { z: 108, s: "Hs", n: "Hassium", cn: "Hassium", type: "metal", state: "solid", mp: "-", bp: "-", p: 7, g: "8B", iupac: 8, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f14 6d6 7s2", noble: "[Rn] 5f14 6d6 7s2" },
    { z: 109, s: "Mt", n: "Meitnerium", cn: "Meitnerium", type: "metal", state: "solid", mp: "-", bp: "-", p: 7, g: "8B", iupac: 9, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f14 6d7 7s2", noble: "[Rn] 5f14 6d7 7s2" },
    { z: 110, s: "Ds", n: "Darmstadtium", cn: "Darmstadtium", type: "metal", state: "solid", mp: "-", bp: "-", p: 7, g: "8B", iupac: 10, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f14 6d8 7s2", noble: "[Rn] 5f14 6d8 7s2" },
    { z: 111, s: "Rg", n: "Roentgenium", cn: "Roentgenium", type: "metal", state: "solid", mp: "-", bp: "-", p: 7, g: "1B", iupac: 11, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f14 6d9 7s2", noble: "[Rn] 5f14 6d9 7s2" },
    { z: 112, s: "Cn", n: "Copernicium", cn: "Copernicium", type: "metal", state: "solid", mp: "-", bp: "-", p: 7, g: "2B", iupac: 12, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f14 6d10 7s2", noble: "[Rn] 5f14 6d10 7s2" },
    { z: 113, s: "Nh", n: "Nihonium", cn: "Nihonium", type: "metal", state: "solid", mp: "-", bp: "-", p: 7, g: "3A", iupac: 13, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f14 6d10 7s2 7p1", noble: "[Rn] 5f14 6d10 7s2 7p1" },
    { z: 114, s: "Fl", n: "Flerovium", cn: "Flerovium", type: "metal", state: "solid", mp: "-", bp: "-", p: 7, g: "4A", iupac: 14, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f14 6d10 7s2 7p2", noble: "[Rn] 5f14 6d10 7s2 7p2" },
    { z: 115, s: "Mc", n: "Moscovium", cn: "Moscovium", type: "metal", state: "solid", mp: "-", bp: "-", p: 7, g: "5A", iupac: 15, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f14 6d10 7s2 7p3", noble: "[Rn] 5f14 6d10 7s2 7p3" },
    { z: 116, s: "Lv", n: "Livermorium", cn: "Livermorium", type: "metal", state: "solid", mp: "-", bp: "-", p: 7, g: "6A", iupac: 16, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f14 6d10 7s2 7p4", noble: "[Rn] 5f14 6d10 7s2 7p4" },
    { z: 117, s: "Ts", n: "Tennessine", cn: "Tennessine", type: "metalloid", state: "solid", mp: "-", bp: "-", p: 7, g: "7A", iupac: 17, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f14 6d10 7s2 7p5", noble: "[Rn] 5f14 6d10 7s2 7p5" },
    { z: 118, s: "Og", n: "Oganesson", cn: "Oganesson", type: "nonmetal", state: "gas", mp: "-", bp: "-", p: 7, g: "8A", iupac: 18, c: "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6 5f14 6d10 7s2 7p6", noble: "[Rn] 5f14 6d10 7s2 7p6" }
];

// ========== [Daily fortune database: final revision] ==========

// ====== Cleanup: remove "模式: 複合" / composite-mode artifacts if present ======
(function cleanupCompositeMode() {
    const candidateNames = [
        "模式: 複合","模式_複合","複合模式",
        "compositeMode","modeComposite","composite"
    ];

    function tryDeleteGlobal(name){
        try{
            if (typeof window !== 'undefined' && window.hasOwnProperty(name)) {
                try { delete window[name]; } catch(e) { window[name] = undefined; }
            }
        }catch(e){}
    }

    function removeDomReferences(){
        const attrSelectors = ['[data-mode]','[data-mode-name]','[data-mode-id]','[data-mode-type]'];
        document.querySelectorAll(attrSelectors.join(',')).forEach(el=>{
            const v = (el.getAttribute('data-mode') || el.getAttribute('data-mode-name') || el.getAttribute('data-mode-id') || el.getAttribute('data-mode-type') || '');
            if (v.indexOf('複合') !== -1 || /composite/i.test(v)) el.remove();
        });
        // remove nodes whose text contains "模式: 複合"
        document.querySelectorAll('body *').forEach(el=>{
            if (el.childElementCount === 0 && /模式[:：]\s*複合/.test((el.textContent||''))) el.remove();
        });
    }

    function removeByNamePatterns(){
        // attempt to clear any globally reachable variables/functions containing keywords
        Object.getOwnPropertyNames(window).forEach(k=>{
            if (/複合|composite/i.test(k)) {
                try { delete window[k]; } catch(e) { window[k] = undefined; }
            }
        });
        // explicit list
        candidateNames.forEach(tryDeleteGlobal);
    }

    function run(){
        try { removeByNamePatterns(); } catch(e){}
        try { if (typeof document !== 'undefined') removeDomReferences(); } catch(e){}
        // As a last resort, try to null any likely variables we can find on window
        try {
            candidateNames.forEach(n=>{
                try {
                    if (window[n] !== undefined) window[n] = undefined;
                } catch(e){}
            });
        } catch(e){}
    }

    if (typeof document !== 'undefined' && document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', run);
    } else {
        // run ASAP
        setTimeout(run, 0);
    }
})();