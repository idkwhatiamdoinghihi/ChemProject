// ==========================================
// data_molecules.js - molecule definitions and helper functions
// ==========================================

let MOLECULE_DB = {};
let MOLECULE_INDEX = {};

// 幾何輔助函式
const di = (e, dist=60) => [{elem:e,x:-dist,y:0,z:0},{elem:e,x:dist,y:0,z:0}];
const getLinear = (c, o, r=70) => [ {elem:c,x:0,y:0,z:0, lpCount:0}, {elem:o,x:-r,y:0,z:0}, {elem:o,x:r,y:0,z:0} ];
const getTrigPlanar = (c, o, r=70) => [ {elem:c,x:0,y:0,z:0, lpCount:0}, {elem:o,x:0,y:r,z:0}, {elem:o,x:r*0.866,y:-r*0.5,z:0}, {elem:o,x:-r*0.866,y:-r*0.5,z:0} ];
const getTetra = (c, o, d=60) => { const r = d / 1.73205; return [ {elem:c,x:0,y:0,z:0, lpCount:0}, {elem:o,x:r,y:-r,z:r}, {elem:o,x:-r,y:r,z:r}, {elem:o,x:-r,y:-r,z:-r}, {elem:o,x:r,y:r,z:-r} ]; };
const getOcta = (c, o, r=65) => [{elem:c,x:0,y:0,z:0, lpCount:0}, {elem:o,x:r,y:0,z:0}, {elem:o,x:-r,y:0,z:0}, {elem:o,x:0,y:r,z:0}, {elem:o,x:0,y:-r,z:0}, {elem:o,x:0,y:0,z:r}, {elem:o,x:0,y:0,z:-r}];
const benzBase=[{x:0,y:70,z:0},{x:60,y:35,z:0},{x:60,y:-35,z:0},{x:0,y:-70,z:0},{x:-60,y:-35,z:0},{x:-60,y:35,z:0}];
function getBenzH(i,s=35){const v=benzBase[i],l=Math.sqrt(v.x**2+v.y**2);return{x:v.x+v.x/l*s,y:v.y+v.y/l*s,z:0};}

// Main registration function (modified: use the passed hybrid string directly)
const addMol = (keysStr, center, hybrid, shape, angle, mp, bp, atoms, bonds, variants = null, desc = null, pg = null) => {
    if (typeof MOLECULE_INDEX === 'undefined') MOLECULE_INDEX = {};
    if (typeof MOLECULE_DB === 'undefined') MOLECULE_DB = {};
    const keys = keysStr.split('|');
    const mainKey = keys[0].trim();
    const mainKeyUpper = mainKey.toUpperCase();
    keys.forEach(k => { MOLECULE_INDEX[k.trim().toUpperCase()] = { key: mainKey, variant: null }; });

    const baseData = { 
        center, 
        hybrid: hybrid, 
        shape: Array.isArray(shape) ? `${shape[0]} (${shape[1]})` : shape, 
        angle, mp, bp, atomsRaw: atoms, bondsRaw: bonds, desc, fullKey: keysStr,
        isMetal: false,
        pg: pg // 確保這裡有接收到傳入的 pg
    };
    
    if (variants) {
        baseData.variants = {};
        for (let vKeyRaw in variants) {
            const uniqueID = vKeyRaw; 
            const vObj = variants[vKeyRaw];
            const vKeys = vKeyRaw.split('|');
            vKeys.forEach(vk => { 
                const cleanKey = vk.trim().toUpperCase();
                if (cleanKey !== mainKeyUpper) { MOLECULE_INDEX[cleanKey] = { key: mainKey, variant: uniqueID }; }
            });
            baseData.variants[uniqueID] = { 
                ...baseData, 
                atomsRaw: vObj.atoms, 
                bondsRaw: vObj.bonds,
                // Core fix: ensure the pg tag can be read from data and passed through
                pg: vObj.pg || baseData.pg || null,
                mp: vObj.mp !== undefined ? vObj.mp : baseData.mp,
                bp: vObj.bp !== undefined ? vObj.bp : baseData.bp,
                desc: vObj.desc !== undefined ? vObj.desc : baseData.desc,
                fullKey: vKeyRaw 
            };
        }
    }
    MOLECULE_DB[mainKey] = baseData;
};

function markReps(atoms, bonds, cnA, elemA, cnB, elemB) {
    const counts = new Array(atoms.length).fill(0);
    bonds.forEach(b => { counts[b[0]]++; counts[b[1]]++; });
    atoms.forEach((a, i) => {
        if (a.elem === elemA && counts[i] === cnA) a.isRepresentative = true;
        else if (a.elem === elemB && counts[i] === cnB) a.isRepresentative = true;
        else a.isRepresentative = false;
    });
}

// NaCl crystal
(function(){
    const sa=[{elem:"Na",x:-40,y:0,z:0,r:20,lpCount:0},{elem:"Cl",x:40,y:0,z:0,r:35,lpCount:0}], sb=[[0,1,"ionic_thin"]];
    const ca=[], cb=[], s=120;
    for(let x=-1;x<=1;x++) for(let y=-1;y<=1;y++) for(let z=-1;z<=1;z++){
        const isNa=(Math.abs(x+y+z)%2!==0);
        ca.push({elem:isNa?"Na":"Cl",x:x*s,y:y*s,z:z*s,r:isNa?18:34,lpCount:0,gx:x,gy:y,gz:z,isRepresentative:(!x&&!y&&!z)});
    }
    for(let i=0;i<ca.length;i++) for(let j=i+1;j<ca.length;j++){
        const dist=Math.abs(ca[i].x-ca[j].x)+Math.abs(ca[i].y-ca[j].y)+Math.abs(ca[i].z-ca[j].z);
        if(Math.abs(dist-s)<1){
            const onFace=(Math.abs(ca[i].gx)===1&&ca[i].gx===ca[j].gx)||(Math.abs(ca[i].gy)===1&&ca[i].gy===ca[j].gy)||(Math.abs(ca[i].gz)===1&&ca[i].gz===ca[j].gz);
            cb.push([i,j,onFace?"ionic_thick":"ionic_thin"]);
        }
    }
    addMol("NaCl|Sodium Chloride|table salt","Na","-","-","-","801","1465",sa,sb,{
        "Simple|Basic unit (ion pair)":{atoms:sa,bonds:sb,hybrid:"-",shape:"-",desc:'<div class="info-section"><div class="info-title">🧂 Overview</div><div class="info-body"><strong>Sodium chloride (NaCl)</strong><br>Commonly known as table salt. In pure form it is a transparent crystal. Widely used as seasoning and preservative.</div></div>'},
        "Crystal|Crystal packing (FCC)":{atoms:ca,bonds:cb,isIonic:true,edgeRelation:"a = 2(r<sub>+</sub> + r<sub>-</sub>)",desc:'<div class="info-section"><div class="info-title">🧊 Crystal properties</div><div class="info-body"><strong>Face-centered cubic (FCC)</strong><br>Sodium chloride has a high melting point (801°C). Each Na+ is surrounded by 6 Cl-, giving a coordination number of 6.<br><span style="color:#facc15">★ Click the center atom to view coordination number.</span></div></div>'}
    },'<div class="info-section"><div class="info-title">🧂 Sodium Chloride</div><div class="info-body">Switch variants to view details.</div></div>');
    if(MOLECULE_DB["NaCl"]?.variants){MOLECULE_DB["NaCl"].variants["Crystal|晶體堆積 (FCC)"].isIonic=true;MOLECULE_DB["NaCl"].variants["Simple|基本單元 (離子對)"].isIonic=true;}
})();

// CsCl crystal
(function(){
    const sa=[{elem:"Cs",x:-45,y:0,z:0,r:26,lpCount:0},{elem:"Cl",x:45,y:0,z:0,r:34,lpCount:0}], sb=[[0,1,"ionic_thin"]];
    const ca=[], cb=[], s=200;
    ca.push({elem:"Cs",x:0,y:0,z:0,r:26,isRepresentative:true});
    [-1,1].forEach(x=>[-1,1].forEach(y=>[-1,1].forEach(z=>{ca.push({elem:"Cl",x:x*s*0.5,y:y*s*0.5,z:z*s*0.5,r:34,isCorner:true}); cb.push([0,ca.length-1,"ionic_thin"]);})));
    for(let i=1;i<ca.length;i++) for(let j=i+1;j<ca.length;j++) if(Math.abs((Math.abs(ca[i].x-ca[j].x)+Math.abs(ca[i].y-ca[j].y)+Math.abs(ca[i].z-ca[j].z))-s)<5) cb.push([i,j,"ionic_thick"]);
    addMol("CsCl|Cesium Chloride","Cs","-","-","-","645","1290",sa,sb,{
        "Simple|Basic unit (ion pair)":{atoms:sa,bonds:sb,hybrid:"-",shape:"-",desc:'<div class="info-section"><div class="info-title">⚛️ Overview</div><div class="info-body"><strong>Cesium chloride (CsCl)</strong><br>Composed of Cs+ and Cl-. Cs+ has a large ionic radius, forming a coordination number of 8.</div></div>'},
        "Crystal|Crystal packing (SC)":{atoms:ca,bonds:cb,isIonic:true,edgeRelation:"√3 a = 2(r⁺+r⁻)",desc:'<div class="info-section"><div class="info-title">🧊 Crystal structure</div><div class="info-body"><strong>Simple cubic (SC)</strong><br>Chloride ions form a simple cubic lattice with cesium ions occupying the body center. Coordination number is 8.<br><span style="color:#facc15">★ Click the central Cs+ to view coordination number.</span></div></div>'}
    },'<div class="info-section"><div class="info-title">🧊 Cesium Chloride</div><div class="info-body">Switch variants to view details.</div></div>');
    if(MOLECULE_DB["CsCl"]?.variants){MOLECULE_DB["CsCl"].variants["Crystal|晶體堆積 (SC)"].isIonic=true;MOLECULE_DB["CsCl"].variants["Simple|基本單元 (離子對)"].isIonic=true;}
})();

// ZnS crystal
(function(){
    const sa=[{elem:"Zn",x:-45,y:0,z:0,r:18,lpCount:0},{elem:"S",x:45,y:0,z:0,r:30,lpCount:0}], sb=[[0,1,"ionic_thin"]];
    const ca=[], cb=[], scale=220, bondDist=scale*0.433; 
    const baseS=[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1],[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5],[0.5,1,0.5],[1,0.5,0.5],[0.5,0.5,1]];
    const baseZn=[[0.25,0.25,0.25],[0.75,0.75,0.25],[0.75,0.25,0.75],[0.25,0.75,0.75]];
    let idx=0;
    baseS.forEach((p,i)=>ca.push({elem:"S",x:(p[0]-0.5)*scale,y:(p[1]-0.5)*scale,z:(p[2]-0.5)*scale,r:28,isCorner:(i<8),idx:idx++}));
    baseZn.forEach(p=>ca.push({elem:"Zn",x:(p[0]-0.5)*scale,y:(p[1]-0.5)*scale,z:(p[2]-0.5)*scale,r:12,isRepresentative:true}));
    for(let i=14;i<ca.length;i++) for(let j=0;j<14;j++) if(Math.abs(Math.sqrt((ca[i].x-ca[j].x)**2+(ca[i].y-ca[j].y)**2+(ca[i].z-ca[j].z)**2)-bondDist)<20) cb.push([i,j,"ionic_thin"]);
    for(let i=0;i<14;i++) for(let j=i+1;j<14;j++) if(ca[i].isCorner&&ca[j].isCorner&&Math.abs(Math.sqrt((ca[i].x-ca[j].x)**2+(ca[i].y-ca[j].y)**2+(ca[i].z-ca[j].z)**2)-scale)<5) cb.push([i,j,"ionic_thick"]);
    addMol("ZnS|Zinc Blende|zinc sulfide","Zn","-","-","-","1185","sublimation",sa,sb,{
        "Simple|Basic unit (ion pair)":{atoms:sa,bonds:sb,hybrid:"-",shape:"-",desc:'<div class="info-section"><div class="info-title">💡 Properties</div><div class="info-body"><strong>Zinc sulfide (ZnS)</strong><br>White to pale yellow powder with luminescent properties. When doped with trace metals it is used in phosphors, glow-in-the-dark paints, and display phosphors.</div></div>'},
        "Crystal|Crystal packing (FCC)":{atoms:ca,bonds:cb,isIonic:true,edgeRelation:"4(r<sub>+</sub> + r<sub>-</sub>) = √3 a",desc:'<div class="info-section"><div class="info-title">💎 Zinc blende (ZnS)</div><div class="info-body">S²⁻ forms a face-centered cubic packing while Zn²⁺ occupies tetrahedral sites.<br><span style="color:#facc15">★ Click any internal Zn to view its coordination number.</span></div></div>'}
    },'<div class="info-section"><div class="info-title">💡 Zinc Sulfide</div><div class="info-body">Switch variants to view details.</div></div>');
    if(MOLECULE_DB["ZnS"]?.variants){MOLECULE_DB["ZnS"].variants["Crystal|晶體堆積 (FCC)"].isIonic=true;MOLECULE_DB["ZnS"].variants["Simple|基本單元 (離子對)"].isIonic=true;}
})();

// CuCl crystal
(function(){
    const sa=[{elem:"Cu",x:-45,y:0,z:0,r:13,lpCount:0},{elem:"Cl",x:45,y:0,z:0,r:27,lpCount:0}], sb=[[0,1,"ionic_thin"]];
    const ca=[], cb=[], scale=220, bondDist=scale*0.433;
    const baseCl=[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1],[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5],[0.5,1,0.5],[1,0.5,0.5],[0.5,0.5,1]];
    const baseCu=[[0.25,0.25,0.25],[0.75,0.75,0.25],[0.75,0.25,0.75],[0.25,0.75,0.75]];
    let clIdx=0;
    baseCl.forEach((p,i)=>ca.push({elem:"Cl",x:(p[0]-0.5)*scale,y:(p[1]-0.5)*scale,z:(p[2]-0.5)*scale,r:27,isCorner:(i<8),idx:clIdx++}));
    baseCu.forEach(p=>ca.push({elem:"Cu",x:(p[0]-0.5)*scale,y:(p[1]-0.5)*scale,z:(p[2]-0.5)*scale,r:13,isRepresentative:true}));
    for(let i=14;i<ca.length;i++) for(let j=0;j<14;j++) if(Math.abs(Math.sqrt((ca[i].x-ca[j].x)**2+(ca[i].y-ca[j].y)**2+(ca[i].z-ca[j].z)**2)-bondDist)<20) cb.push([i,j,"ionic_thin"]);
    for(let i=0;i<14;i++) for(let j=i+1;j<14;j++) if(ca[i].isCorner&&ca[j].isCorner&&Math.abs(Math.sqrt((ca[i].x-ca[j].x)**2+(ca[i].y-ca[j].y)**2+(ca[i].z-ca[j].z)**2)-scale)<5) cb.push([i,j,"ionic_thick"]);
    addMol("CuCl|Copper(I) chloride|Nantokite","Cu","-","-","-","430","1490",sa,sb,{
        "Simple|Basic unit (ion pair)":{atoms:sa,bonds:sb,hybrid:"-",shape:"-",desc:'<div class="info-section"><div class="info-title">🔸 Overview</div><div class="info-body"><strong>Copper(I) chloride (CuCl)</strong><br>White solid, sparingly soluble in water. Structure similar to zinc blende (ZnS).</div></div>'},
        "Crystal|Crystal packing (FCC)":{atoms:ca,bonds:cb,isIonic:true,edgeRelation:"4(r<sub>+</sub> + r<sub>-</sub>) = √3 a",desc:'<div class="info-section"><div class="info-title">🧊 Crystal structure</div><div class="info-body"><strong>Face-centered cubic (FCC)</strong><br>Structure is analogous to zinc blende; Cl- packs in an FCC lattice and Cu+ occupies tetrahedral holes.<br><span style="color:#facc15">★ Click any internal Cu⁺ to view coordination number.</span></div></div>'}
    },'<div class="info-section"><div class="info-title">🔸 Copper(I) chloride</div><div class="info-body">Switch variants to view details.</div></div>');
    if(MOLECULE_DB["CuCl"]?.variants){MOLECULE_DB["CuCl"].variants["Crystal|晶體堆積 (FCC)"].isIonic=true;MOLECULE_DB["CuCl"].variants["Simple|基本單元 (離子對)"].isIonic=true;}
})();

// TiO2 crystal
(function(){
    const sa=[{elem:"Ti",x:0,y:0,z:0,r:11,lpCount:0},{elem:"O",x:50,y:0,z:0,r:21,lpCount:0},{elem:"O",x:-50,y:0,z:0,r:21,lpCount:0}], sb=[[0,1,"ionic_thin"],[0,2,"ionic_thin"]];
    const ca=[], cb=[], scale=180, c_ratio=0.65, u=0.3;
    const baseTi=[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1],[0.5,0.5,0.5]];
    const baseO=[[u,u,0],[1-u,1-u,0],[u,u,1],[1-u,1-u,1],[0.5+u,0.5-u,0.5],[0.5-u,0.5+u,0.5]];
    let tiIdx=0;
    baseTi.forEach((p,i)=>ca.push({elem:"Ti",x:(p[0]-0.5)*scale,y:(p[1]-0.5)*scale,z:(p[2]-0.5)*scale*c_ratio,r:11,isCorner:(i<8),idx:tiIdx++,isRepresentative:(i===8)}));
    baseO.forEach(p=>ca.push({elem:"O",x:(p[0]-0.5)*scale,y:(p[1]-0.5)*scale,z:(p[2]-0.5)*scale*c_ratio,r:21,isCorner:false}));
    for(let i=0;i<ca.length;i++) for(let j=i+1;j<ca.length;j++){
        if(ca[i].elem===ca[j].elem) continue;
        if(Math.sqrt((ca[i].x-ca[j].x)**2+(ca[i].y-ca[j].y)**2+(ca[i].z-ca[j].z)**2)<scale*0.75) cb.push([i,j,"ionic_thin"]);
    }
    for(let i=0;i<8;i++) for(let j=i+1;j<8;j++){
        const dx=Math.abs(ca[i].x-ca[j].x), dy=Math.abs(ca[i].y-ca[j].y), dz=Math.abs(ca[i].z-ca[j].z);
        if((Math.abs(dx-scale)<5&&dy<5&&dz<5)||(Math.abs(dy-scale)<5&&dx<5&&dz<5)||(Math.abs(dz-scale*c_ratio)<5&&dx<5&&dy<5)) cb.push([i,j,"ionic_thick"]);
    }
    addMol("TiO2|Rutile|Titanium dioxide","Ti","-","-","-","1843","2972",sa,sb,{
        "Simple|Basic unit":{atoms:sa,bonds:sb,hybrid:"-",shape:"-",desc:'<div class="info-section"><div class="info-title">⬜ Overview</div><div class="info-body"><strong>Titanium dioxide (TiO₂)</strong><br>White powder widely used as a white pigment, in sunscreens, and as a photocatalyst.</div></div>'},
        "Crystal|Crystal packing (Tetragonal)":{atoms:ca,bonds:cb,isIonic:true,edgeRelation:"complex geometry",desc:'<div class="info-section"><div class="info-title">🧊 Crystal structure</div><div class="info-body"><strong>Tetragonal (Rutile)</strong><br>Titanium ions occupy body-centered and corner positions while oxygen sits on faces. Ti⁴⁺ has coordination number 6 (octahedral); O²⁻ has CN=3 (planar triangle).<br><span style="color:#facc15">★ Click the body-centered Ti⁴⁺ to view its coordination number.</span></div></div>'}
    },'<div class="info-section"><div class="info-title">⬜ Rutile</div><div class="info-body">Switch variants to view details.</div></div>');
    if(MOLECULE_DB["TiO2"]?.variants){MOLECULE_DB["TiO2"].variants["Crystal|晶體堆積 (Tetragonal)"].isIonic=true;MOLECULE_DB["TiO2"].variants["Simple|基本單元"].isIonic=true;}
})();

// Cu2O crystal
(function(){
    const scale=180, baseO=[[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1],[0.5,0.5,0.5]], baseCu=[[0.25,0.25,0.25],[0.75,0.75,0.25],[0.75,0.25,0.75],[0.25,0.75,0.75]];
    const ca=[...baseO.map((p,i)=>({elem:"O",x:(p[0]-0.5)*scale,y:(p[1]-0.5)*scale,z:(p[2]-0.5)*scale,r:21,isCorner:i<8,isRepresentative:i===8})),...baseCu.map(p=>({elem:"Cu",x:(p[0]-0.5)*scale,y:(p[1]-0.5)*scale,z:(p[2]-0.5)*scale,r:13,isRepresentative:true}))];
    const cb=[];
    for(let i=0;i<ca.length;i++) for(let j=i+1;j<ca.length;j++){
        const d=Math.hypot(ca[i].x-ca[j].x,ca[i].y-ca[j].y,ca[i].z-ca[j].z);
        if(ca[i].elem!==ca[j].elem&&Math.abs(d-scale*0.433)<20) cb.push([i,j,"ionic_thin"]);
        if(ca[i].isCorner&&ca[j].isCorner&&Math.abs(d-scale)<5) cb.push([i,j,"ionic_thick"]);
    }
    const sa=[{elem:"O",x:0,y:0,z:0,r:21},{elem:"Cu",x:50,y:0,z:0,r:13},{elem:"Cu",x:-50,y:0,z:0,r:13}], sb=[[0,1,"ionic_thin"],[0,2,"ionic_thin"]];
    addMol("Cu2O|Cuprite|Copper(I) oxide","Cu","-","-","-","1235","1800",sa,sb,{
        "Simple|Basic unit":{atoms:sa,bonds:sb,hybrid:"-",shape:"-",desc:'<div class="info-section"><div class="info-title">🔴 Overview</div><div class="info-body"><strong>Copper(I) oxide (Cu₂O)</strong><br>Red solid. Cu⁺ shows linear coordination (CN=2), while O²⁻ has tetrahedral coordination (CN=4).</div></div>'},
        "Crystal|Crystal packing (Cubic)":{atoms:ca,bonds:cb,isIonic:true,edgeRelation:"complex geometry",desc:'<div class="info-section"><div class="info-title">🧊 Crystal structure</div><div class="info-body"><strong>Cuprite structure</strong><br>Oxygen ions (red) form a body-centered cubic arrangement; copper ions (orange) sit at midpoints between oxygens.<br>• Click a <strong>red oxygen</strong> (body center) to see CN=4.<br>• Click any <strong>orange copper</strong> to see CN=2.</div></div>'}
    },'<div class="info-section"><div class="info-title">🔴 Cuprite</div><div class="info-body">Switch variants to view details.</div></div>');
    if(MOLECULE_DB["Cu2O"]?.variants){MOLECULE_DB["Cu2O"].variants["Crystal|晶體堆積 (Cubic)"].isIonic=true;MOLECULE_DB["Cu2O"].variants["Simple|基本單元"].isIonic=true;}
})();











// ==========================================
// 金屬晶體生成模組 (CN=12 延伸增強版)
// ==========================================

function ensureElement(elem, defaultColor, defaultR) {
    if (typeof ELEMENT_PROPS !== 'undefined' && !ELEMENT_PROPS[elem]) {
        ELEMENT_PROPS[elem] = { ve: 1, c3d: defaultColor, r3d: defaultR, lp: 0, mass: 0, en: 0 };
    }
}

// 1. 簡單立方 (SC) 
function addMetal_SC(elem, name, mp, bp, scale=160) {
    ensureElement(elem, "#ab5c00", 28);
    const atoms = []; const bonds = [];
    for (let x = 0; x <= 1; x++) {
        for (let y = 0; y <= 1; y++) {
            for (let z = 0; z <= 1; z++) {
                atoms.push({ elem: elem, x: (x-0.5)*scale, y: (y-0.5)*scale, z: (z-0.5)*scale, r: 28, isRepresentative: true });
            }
        }
    }
    for (let i = 0; i < atoms.length; i++) {
        for (let j = i + 1; j < atoms.length; j++) {
            const d = Math.sqrt((atoms[i].x-atoms[j].x)**2 + (atoms[i].y-atoms[j].y)**2 + (atoms[i].z-atoms[j].z)**2);
            if (Math.abs(d - scale) < 10) bonds.push([i, j, "ionic_thick"]);
        }
    }
    addMol(`${elem}|${name}`, "Metal", "簡單立方堆積 (SC)", "52.4%", "6", mp, bp, atoms, bonds, null,
        `<div class="info-section"><div class="info-title">📦 簡單立方 (SC)</div><div class="info-body">金屬範例：<strong>${elem}</strong>。<br>空間利用率 52.4%。原子僅位於立方體頂點，沿著邊長互相接觸。</div></div>`);
    if(MOLECULE_DB[elem]) { MOLECULE_DB[elem].isIonic = true; MOLECULE_DB[elem].isMetal = true; MOLECULE_DB[elem].edgeRelation = "a = 2r"; }
}

// 2. 體心立方 (BCC)
function addMetal_BCC(elem, name, mp, bp, scale=200) {
    ensureElement(elem, "#9ca3af", 24); 
    const atoms = []; const bonds = [];
    const h = scale / 2;
    atoms.push({ elem: elem, x: 0, y: 0, z: 0, r: 24, isRepresentative: true });
    const pts = [[-h,-h,-h],[h,-h,-h],[h,h,-h],[-h,h,-h],[-h,-h,h],[h,-h,h],[h,h,h],[-h,h,h]];
    pts.forEach(p => atoms.push({ elem: elem, x: p[0], y: p[1], z: p[2], r: 24 }));
    for(let i=1; i<=8; i++) bonds.push([0, i, "ionic_thin"]);
    const cubeEdges = [[1,2],[2,3],[3,4],[4,1],[5,6],[6,7],[7,8],[8,5],[1,5],[2,6],[3,7],[4,8]];
    cubeEdges.forEach(e => bonds.push([e[0], e[1], "ionic_thick"]));
    addMol(`${elem}|${name}`, "Metal", "體心立方堆積 (BCC)", "68%", "8", mp, bp, atoms, bonds, null,
        `<div class="info-section"><div class="info-title">🧊 體心立方 (BCC)</div><div class="info-body">金屬範例：<strong>${elem}</strong>。<br>空間利用率 68%。原子位於角落與體中心，沿著體對角線互相接觸。</div></div>`);
    if(MOLECULE_DB[elem]) { MOLECULE_DB[elem].isIonic = true; MOLECULE_DB[elem].isMetal = true; MOLECULE_DB[elem].edgeRelation = "√3 a = 4r"; }
}

// 3. 面心立方 (FCC) - 升級為 5-4-5-4 堆積 (展示 CN=12)
function addMetal_FCC(elem, name, mp, bp, scale=200) {
    ensureElement(elem, "#d1d5db", 22);
    const atoms = []; const bonds = [];
    const h = scale / 2;
    
    // 定義四層：L1(5) -> L2(4) -> L3(5) -> L4(4)
    // 我們將座標中心設在 L3 的中心原子，方便旋轉觀察
    const addLayer5 = (z, isMain) => {
        const s = atoms.length;
        // 中心
        atoms.push({ elem: elem, x: 0, y: 0, z: z, r: 22, isRepresentative: isMain });
        // 四個角
        const corners = [[-h,-h,z],[h,-h,z],[h,h,z],[-h,h,z]];
        corners.forEach(p => atoms.push({ elem: elem, x: p[0], y: p[1], z: p[2], r: 22 }));
        // 只有主要晶胞層 (L1到L3) 有框
        if (isMain || z < scale) {
            bonds.push([s+1, s+2, "ionic_thick"], [s+2, s+3, "ionic_thick"], [s+3, s+4, "ionic_thick"], [s+4, s+1, "ionic_thick"]);
        }
    };

    const addLayer4 = (z, isExtended) => {
        const s = atoms.length;
        // 四個面心
        const faces = [[0,-h,z],[h,0,z],[0,h,z],[-h,0,z]];
        faces.forEach(p => atoms.push({ elem: elem, x: p[0], y: p[1], z: p[2], r: 22 }));
        // 垂直柱子 (只連接 L1-L3 核心)
        if (!isExtended) {
            // 此處邏輯由後續接觸線處理
        }
    };

    addLayer5(-scale, false); // L1 (底部)
    addLayer4(-h, false);     // L2
    addLayer5(0, true);       // L3 (核心層，設為座標 0)
    addLayer4(h, true);       // L4 (延伸層)

    // 建立所有原子間的接觸線 (距離為 0.707a)
    const contactDist = scale * 0.707;
    for (let i = 0; i < atoms.length; i++) {
        for (let j = i + 1; j < atoms.length; j++) {
            const d = Math.sqrt((atoms[i].x-atoms[j].x)**2 + (atoms[i].y-atoms[j].y)**2 + (atoms[i].z-atoms[j].z)**2);
            if (Math.abs(d - contactDist) < 10) {
                // 判斷是否屬於延伸層 (L4) 的連線
                const isExt = (atoms[i].z > 5 || atoms[j].z > 5);
                bonds.push([i, j, isExt ? "ionic_thin" : "ionic_thin"]); 
            }
        }
    }

    // 建立核心晶胞的垂直粗框線 (L1 到 L3)
    const coreCorners = [[1,10],[2,11],[3,12],[4,13]]; // L1 到 L3 的頂點對應
    coreCorners.forEach(e => bonds.push([e[0], e[1], "ionic_thick"]));

    addMol(`${elem}|${name}`, "Metal", "面心立方堆積 (FCC)", "74%", "12", mp, bp, atoms, bonds, null,
        `<div class="info-section"><div class="info-title">✨ 面心立方 (FCC)</div><div class="info-body">金屬範例：<strong>${elem}</strong>。<br>空間利用率 74%。模型展示了 5-4-5-4 的四層堆積。<br><span style="color:#facc15">★ 點擊第三層中心原子，可見其配位數為 12 (同層4, 下層4, 上層4)。</span></div></div>`);
    if(MOLECULE_DB[elem]) { MOLECULE_DB[elem].isIonic = true; MOLECULE_DB[elem].isMetal = true; MOLECULE_DB[elem].edgeRelation = "√2 a = 4r"; }
}

// 4. 六方最密堆積 (HCP) - 修正比例與配位數版
function addMetal_HCP(elem, name, mp, bp, scale=140) {
    ensureElement(elem, "#e5e7eb", 22);
    const atoms = []; const bonds = [];
    
    // a = scale (原子間距，即底面六角形的邊長)
    // h = 層與層之間的垂直距離 (理想比例為 sqrt(2/3) * a ≈ 0.8165a)
    const h = scale * 0.8165; 
    const r = 22;

    // A 層生成器 (中心 + 6 顆環繞)
    const getLayerA = (z) => [
        {x:0, y:0, z:z}, // 中心
        {x:scale, y:0, z:z}, {x:scale*0.5, y:scale*0.866, z:z}, {x:-scale*0.5, y:scale*0.866, z:z},
        {x:-scale, y:0, z:z}, {x:-scale*0.5, y:-scale*0.866, z:z}, {x:scale*0.5, y:-scale*0.866, z:z}
    ];

    // B 層生成器 (填入 A 層空隙的 3 顆)
    const getLayerB = (z) => [
        {x:scale*0.5, y:scale*0.288, z:z}, 
        {x:-scale*0.5, y:scale*0.288, z:z}, 
        {x:0, y:-scale*0.577, z:z}
    ];

    // --- 建立四層堆積 A1-B1-A2-B2 ---
    // 為了讓中心原子在座標原點，我們這樣對齊：
    const l1 = getLayerA(-2 * h);   // Index 0-6 (A1 最底層)
    const l2 = getLayerB(-h);       // Index 7-9 (B1)
    const l3 = getLayerA(0);        // Index 10-16 (A2 核心主角層)
    const l4 = getLayerB(h);        // Index 17-19 (B2 延伸層)

    [...l1, ...l2, ...l3, ...l4].forEach((p, i) => {
        atoms.push({
            elem: elem, ...p, r: r, lpCount: 0,
            // 將第三層的中心原子 (Index 10) 設為主角
            isRepresentative: (i === 10) 
        });
    });

    // --- 建立鍵結邏輯 ---
    for (let i = 0; i < atoms.length; i++) {
        for (let j = i + 1; j < atoms.length; j++) {
            const d = Math.sqrt((atoms[i].x-atoms[j].x)**2 + (atoms[i].y-atoms[j].y)**2 + (atoms[i].z-atoms[j].z)**2);
            
            // 距離約等於 scale (a) 的判定為鄰居
            if (d > 10 && d < scale * 1.1) {
                // 判斷是否為延伸層 (L4 / Index 17-19) 的連線
                const isExt = (atoms[i].z > h/2 || atoms[j].z > h/2);
                
                // 1. 同層內部的連線 (六角外框與內部輻射)
                if (Math.abs(atoms[i].z - atoms[j].z) < 1) {
                    const isCenter = (Math.abs(atoms[i].x) < 1 && Math.abs(atoms[i].y) < 1) || 
                                     (Math.abs(atoms[j].x) < 1 && Math.abs(atoms[j].y) < 1);
                    
                    if (isCenter) {
                        bonds.push([i, j, "ionic_thin"]); // 內部輻射用細線
                    } else {
                        // 外部六角框：L1到L3用粗框，延伸層L4用細線
                        bonds.push([i, j, isExt ? "ionic_thin" : "ionic_thick"]);
                    }
                } 
                // 2. 層與層之間的連線 (CN=12 的斜向接觸)
                else {
                    bonds.push([i, j, "ionic_thin"]);
                }
            }
        }
    }

    // --- 核心六角柱的「垂直」稜線 (L1 頂點對應到 L3 頂點) ---
    // 讓主要的晶胞框架看起來像一個完整的六角柱
    for (let i = 1; i <= 6; i++) {
        bonds.push([i, i + 10, "ionic_thick"]);
    }

    addMol(`${elem}|${name}`, "Metal", "六方最密堆積 (HCP)", "74%", "12", mp, bp, atoms, bonds, null,
        `<div class="info-section"><div class="info-title">🛑 六方最密堆積 (HCP)</div><div class="info-body">金屬範例：<strong>${elem} (如鎂、鋅)</strong>。<br>利用率 74%。模型展示 A-B-A-B 四層堆積，延伸出一層三角形 B 層。<br><span style="color:#facc15">★ 點擊第三層中心原子，可見配位數為 12 (同層6，下層3，上層3)。</span></div></div>`);
    
    if(MOLECULE_DB[elem]) { 
        MOLECULE_DB[elem].isIonic = true; 
        MOLECULE_DB[elem].isMetal = true; 
        MOLECULE_DB[elem].edgeRelation = "c ≈ 1.633 a";
    }
}

// 執行金屬生成
// 1A 族 (BCC)
addMetal_BCC("Li", "鋰", "180.5", "1342");
addMetal_BCC("Na", "鈉", "97.8", "883");
addMetal_BCC("K",  "鉀", "63.5", "759");
addMetal_BCC("Rb", "銣", "39.3", "688");
addMetal_BCC("Cs", "銫", "28.4", "671");

// 2A 族
addMetal_HCP("Be", "鈹", "1287", "2469"); // HCP
addMetal_HCP("Mg", "鎂", "650", "1090");  // HCP
addMetal_FCC("Ca", "鈣", "842", "1484");  // FCC
addMetal_FCC("Sr", "鍶", "777", "1382");  // FCC
addMetal_BCC("Ba", "鋇", "727", "1897");  // BCC

// 其他
addMetal_SC("Po", "釙", "254", "962");
addMetal_BCC("Fe", "鐵 (α)", "1538", "2861");
addMetal_FCC("Cu", "銅", "1085", "2562");
addMetal_FCC("Ag", "銀", "961.8", "2162");
addMetal_FCC("Au", "金", "1064", "2970");
addMetal_FCC("Al", "鋁", "660", "2519");
addMetal_HCP("Zn", "鋅", "419.5", "907");
addMetal_HCP("Ti", "鈦", "1668", "3287");






/*
 ==========================================================================
 ★ 視覺鍵長標準參考表 (Visual Bond Length Standards) v16.2
 ==========================================================================
 基準：以 1,2-二氯丙烷為錨點 (C-C ~ 70, C-H ~ 50, C-Cl ~ 75)
 
 [1] 原子視覺半徑貢獻 (Base Radius Contribution)
 --------------------------------------------------------------------------
  - H (氫) .................... 15  (最小，確保緊湊)
  - Row 2 (C, N, O, F) ........ 35  (基準)
  - Row 3 (Si, P, S, Cl) ...... 40  (略大)
  - Row 4 (Br) ................ 45
  - Row 5 (I, Xe) ............. 50  (最大)

 [2] 鍵級修正係數 (Bond Order Multiplier)
 --------------------------------------------------------------------------
  - 單鍵 (Single) ............. x 1.00
  - 雙鍵 (Double) ............. x 0.90
  - 參鍵 (Triple) ............. x 0.85

 [3] 常見鍵長計算範例 (Calculated Examples)
 --------------------------------------------------------------------------
  Type      Calc (R1 + R2) * Multiplier      Final Value
  -------   ---------------------------      -----------
  H-H       (15 + 15) * 1.0                  30
  C-H       (35 + 15) * 1.0                  50  (基準)
  N-H       (35 + 15) * 1.0                  50
  O-H       (35 + 15) * 1.0                  50
  P-H       (40 + 15) * 1.0                  55

  C-C       (35 + 35) * 1.0                  70  (基準)
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






// ==========================================
// [整理後] 資料注入區 (v14.0 含熔沸點數據)
// ==========================================

/// --- 1. 基礎元素與雙原子分子 (鍵長修正: H=15, 2nd=35, 3rd=40, 4th=45, 5th=50 | Double x0.9, Triple x0.85) ---
const diatomicNames = {'H': '氫|氫氣', 'N': '氮|氮氣', 'O': '氧|氧氣', 'F': '氟|氟氣', 'Cl': '氯|氯氣', 'Br': '溴', 'I': '碘'};
const diatomicProps = {'H': {mp: "-259.2", bp: "-252.9"}, 'N': {mp: "-210.0", bp: "-195.8"}, 'O': {mp: "-218.8", bp: "-183.0"}, 'F': {mp: "-219.7", bp: "-188.1"}, 'Cl': {mp: "-101.5", bp: "-34.0"}, 'Br': {mp: "-7.2", bp: "58.8"}, 'I': {mp: "113.7", bp: "184.3"}};
addMol("H2|氫氣|氫", "雙原子", "s-s", ["直線型", "Linear"], "-", "-259.2", "-252.9", [{elem:"H",x:-15,y:0,z:0},{elem:"H",x:15,y:0,z:0}], [[0,1,"single"]]);
addMol("N2|氮氣|氮", "雙原子", "sp", ["直線型", "Linear"], "-", "-210.0", "-195.8", [{elem:"N",x:-30,y:0,z:0},{elem:"N",x:30,y:0,z:0}], [[0,1,"triple"]]);
addMol("O2|氧氣|氧", "雙原子", "sp²", ["直線型", "Linear"], "-", "-218.8", "-183.0", [{elem:"O",x:-32,y:0,z:0},{elem:"O",x:32,y:0,z:0}], [[0,1,"double"]]);
addMol("F2|氟氣|氟", "雙原子", "sp³", ["直線型", "Linear"], "-", "-219.7", "-188.1", [{elem:"F",x:-35,y:0,z:0},{elem:"F",x:35,y:0,z:0}], [[0,1,"single"]]);
addMol("Cl2|氯氣|氯", "雙原子", "sp³", ["直線型", "Linear"], "-", "-101.5", "-34.0", [{elem:"Cl",x:-40,y:0,z:0},{elem:"Cl",x:40,y:0,z:0}], [[0,1,"single"]]);
addMol("Br2|溴", "雙原子", "sp³", ["直線型", "Linear"], "-", "-7.2", "58.8", [{elem:"Br",x:-45,y:0,z:0},{elem:"Br",x:45,y:0,z:0}], [[0,1,"single"]]);
addMol("I2|碘", "雙原子", "sp³", ["直線型", "Linear"], "-", "113.7", "184.3", [{elem:"I",x:-50,y:0,z:0},{elem:"I",x:50,y:0,z:0}], [[0,1,"single"]]);
["H2", "N2", "O2", "F2", "Cl2", "Br2", "I2"].forEach(key => {
    if (MOLECULE_DB[key]) MOLECULE_DB[key].pg = "Dinfh";
});

addMol("CO|一氧化碳", "雙原子", "sp", ["直線型","Linear"], "-", "-205.0", "-191.5", [{elem:"C",x:-30,y:0,z:0,lp3d:[{x:-1,y:0,z:0}]}, {elem:"O",x:33,y:0,z:0,lp3d:[{x:1,y:0,z:0}]}], [[1,0,"coordinate_triple"]]);
addMol("NO|一氧化氮", "雙原子", "sp²", ["直線型","Linear"], "-", "-164", "-152", [{elem:"N",x:-32,y:0,z:0,radical:true,lp3d:[{x:-1.2,y:1.0,z:0.35},{x:-1.2,y:1.0,z:-0.35},{x:-1.2,y:-1.0,z:0}]},{elem:"O",x:32,y:0,z:0}], [[0,1,"double"]]);
addMol("CN-|氰根|氰離子", "雙原子", "sp", ["直線型","Linear"], "-", "-", "-", [{elem:"C",x:-30,y:0,z:0,lp3d:[{x:-1,y:0,z:0}]},{elem:"N",x:30,y:0,z:0,lp3d:[{x:1,y:0,z:0}]}], [[0,1,"triple"]]);
addMol("O22-|過氧根離子", "O", "sp³", ["直線型","Linear"], "180°", "-", "-", [{elem:"O",x:-35,y:0,z:0,lp3d:[{x:-1,y:1.5,z:0},{x:-1,y:-0.75,z:1.3},{x:-1,y:-0.75,z:-1.3}]}, {elem:"O",x:35,y:0,z:0,lp3d:[{x:1,y:1.5,z:0},{x:1,y:-0.75,z:1.3},{x:1,y:-0.75,z:-1.3}]}], [[0,1]]);
addMol("C22-|碳化物離子", "C", "sp", ["直線型","Linear"], "180°", "-", "-", [{elem:"C",x:-30,y:0,z:0,lp3d:[{x:-1,y:0,z:0}]}, {elem:"C",x:30,y:0,z:0,lp3d:[{x:1,y:0,z:0}]}], [[0,1,"triple"]]);

// --- 2. 鹵化氫 (HX) ---
addMol("HF|氟化氫", "雙原子", "sp³", ["直線型","Linear"], "-", "-83.6", "19.5", [{elem:"F",x:-25,y:0,z:0}, {elem:"H",x:25,y:0,z:0}], [[0,1]]);
addMol("HCl|氯化氫", "雙原子", "sp³", ["直線型","Linear"], "-", "-114.2", "-85.1", [{elem:"Cl",x:-28,y:0,z:0}, {elem:"H",x:28,y:0,z:0}], [[0,1]]);
addMol("HBr|溴化氫", "雙原子", "sp³", ["直線型","Linear"], "-", "-86.8", "-66.4", [{elem:"Br",x:-30,y:0,z:0}, {elem:"H",x:30,y:0,z:0}], [[0,1]]);
addMol("HI|碘化氫", "雙原子", "sp³", ["直線型","Linear"], "-", "-50.8", "-35.4", [{elem:"I",x:-33,y:0,z:0}, {elem:"H",x:33,y:0,z:0}], [[0,1]]);

// 批次設定異核雙原子分子/離子為 Cinfv (直線非對稱)
["CO", "NO", "CN-", "HF", "HCl", "HBr", "HI"].forEach(key => {
    if (MOLECULE_DB[key]) MOLECULE_DB[key].pg = "Cinfv";
});

// 批次設定同核雙原子離子為 Dinfh (直線中心對稱)
["O22-", "C22-"].forEach(key => {
    if (MOLECULE_DB[key]) MOLECULE_DB[key].pg = "Dinfh";
});

// --- 3. 常見無機分子 (H2O, NH3, CH4 等) ---
addMol("CH4|甲烷", "C", "sp³", ["四面體","Tetrahedral"], "109.5°", "-182.5", "-161.5", getTetra("C","H", 50), [[0,1],[0,2],[0,3],[0,4]], null, 
    `<div class="info-section">
        <div class="info-title">⚗️ 物質性質</div>
        <div class="info-body">
            <span class="highlight-title">1. 立體結構：</span>中心碳原子採取 <strong>sp³ 混成軌域</strong>。由於周圍連接四個相同的氫原子且無孤對電子，四個 C-H 鍵之間的電子斥力完全均等，構成了完美的<strong>正四面體</strong>結構，鍵角為 <strong>109.5°</strong>。<br>
            <span class="highlight-title">2. 物理性質：</span>常溫常壓下為無色、無味、無毒的氣體（家用天然氣的臭味是為了安全而添加的硫醇）。屬於完全對稱的<strong>非極性分子</strong>，難溶於水。由於分子量小且分子間僅有微弱的<strong>凡得瓦力</strong>（倫敦分散力），因此熔沸點極低。<br>
            <span class="highlight-title">3. 化學性質：</span>化學性質相當穩定，在一般條件下不與強酸、強鹼或強氧化劑反應。具有可燃性，在空氣中完全燃燒生成二氧化碳與水；在紫外線光照下，可與鹵素（如氯氣）發生連鎖的<strong>自由基取代反應</strong>。
        </div>
    </div>
    <div class="info-section" style="margin-top: 12px; border-top: 1px dashed rgba(255,255,255,0.2); padding-top: 10px;">
        <div class="info-title">🏭 生活應用</div>
        <div class="info-body">
            <span class="highlight-title">1. 潔淨能源 (天然氣)：</span>甲烷是<strong>天然氣</strong>的主要成分 (含量約 90% 以上)。其氫碳比 (H/C ratio) 是所有烴類中最高的，因此燃燒時產生的單位熱值極高，且碳排放量遠低於煤炭與石油，是現代發電與家庭烹飪的重要燃料。<br>
            <span class="highlight-title">2. 未來能源 (可燃冰)：</span>在深海高壓低溫的環境下，甲烷分子會被水分子包覆，形成籠狀結晶結構的<strong>「甲烷水合物」</strong>。外觀晶瑩剔透像冰塊，卻可以直接點火燃燒，其蘊藏量極大，被視為未來最具潛力的戰略能源。<br>
            <span class="highlight-title">3. 溫室氣體效應：</span>雖然大氣中含量遠低於二氧化碳，但甲烷的<strong>全球暖化潛勢 (GWP)</strong> 約是 CO₂ 的 25 倍。這是因為其分子結構中 C-H 鍵的特定震動模式，能非常有效地吸收地表反射的紅外線輻射熱，是造成氣候變遷的關鍵氣體之一。
        </div>
    </div>`, "Td");

addMol("SiH4|矽烷", "Si", "sp³", ["四面體","Tetrahedral"], "109.5°", "-185", "-112", getTetra("Si","H", 55), [[0,1],[0,2],[0,3],[0,4]], null, null, "Td");
addMol("NH3|氨系列", "N", "sp³", ["角錐形","Pyramidal"], "106.7°", "-77.7", "-33.3", [], [], { "NH3|氨|氨氣": {pg: "C3v", mp: "-77.7", bp: "-33.3", desc: "<strong>氨 (Ammonia)</strong><br>三角錐形，具有一對孤對電子，為弱鹼。", atoms: [{elem:"N",x:0,y:10,z:0,lpCount:1}, {elem:"H",x:0,y:-25,z:40}, {elem:"H",x:35,y:-25,z:-20}, {elem:"H",x:-35,y:-25,z:-20}], bonds: [[0,1],[0,2],[0,3]] }, "NH4+|銨根離子|銨離子|銨根": {pg: "Td", mp: "-", bp: "-", desc: "<strong>銨離子</strong><br>正四面體結構，是氨氣與氫離子結合的產物。", atoms: getTetra("N","H", 50), bonds: [[0,1],[0,2],[0,3],[0,4]] }, "NH2-|胺基陰離子|胺基負離子": {pg: "C2v", mp: "-", bp: "-", desc: "<strong>胺基負離子</strong><br>氨失去一個質子後的強鹼性陰離子，V型結構，有兩對孤對電子。", atoms: [{elem:"N",x:0,y:5,z:0,lpCount:2},{elem:"H",x:35,y:-30,z:0},{elem:"H",x:-35,y:-30,z:0}], bonds: [[0,1],[0,2]] }});
addMol("PH3|磷化氫系列", "P", "sp³", ["角錐形","Pyramidal"], "93.3°", "-133.8", "-87.7", [], [], { "PH3|磷化氫": {pg: "C3v", mp: "-133.8", bp: "-87.7", desc: "<strong>磷化氫</strong><br>劇毒氣體，鍵角接近90度(p軌域特性)，但VSEPR視為sp³。", atoms: [{elem:"P",x:0,y:15,z:0,lpCount:1}, {elem:"H",x:0,y:-30,z:45}, {elem:"H",x:39,y:-30,z:-22}, {elem:"H",x:-39,y:-30,z:-22}], bonds: [[0,1],[0,2],[0,3]] }, "PH4+|鏻離子": {pg: "Td", mp: "-", bp: "-", desc: "<strong>鏻離子</strong><br>結構類似銨根，由膦與氫離子形成。", atoms: getTetra("P","H", 55), bonds: [[0,1],[0,2],[0,3],[0,4]] }});
// AA
addMol("H2O|水系列", "O", "sp³", ["角形","Bent"], "104.5°", "0.0", "100.0", [], [], { "H2O|水|水分子": {pg: "C2v",mp: "0.0", bp: "100.0", desc: "<strong>水</strong><br>生命的基石，V型結構，中心氧原子有兩對孤對電子。", atoms: [{elem:"O",x:0,y:5,z:0,lpCount:2}, {elem:"H",x:38,y:-28,z:0}, {elem:"H",x:-38,y:-28,z:0}], bonds: [[0,1],[0,2]] }, "H3O+|水合氫離子|鋞離子": { pg: "C3v",mp: "-", bp: "-", desc: "<strong>水合氫離子</strong><br>水中氫離子的實際存在形式，三角錐形。", atoms: [{elem:"O",x:0,y:10,z:0,lpCount:1}, {elem:"H",x:0,y:-25,z:40}, {elem:"H",x:35,y:-25,z:-20}, {elem:"H",x:-35,y:-25,z:-20}], bonds: [[0,1],[0,2],[0,3]] }, "OH-|氫氧根|氫氧根離子": { pg: "Cinfv",mp: "-", bp: "-", desc: "<strong>氫氧根</strong><br>強鹼的特徵離子，氧原子周圍有三對孤對電子，帶負電。", atoms: [{elem:"O",x:-20,y:0,z:0,lpCount:3},{elem:"H",x:25,y:0,z:0}], bonds: [[0,1]] }});


addMol("H2S|硫化氫系列", "S", "sp³", ["角形","Bent"], "92.1°", "-85.5", "-60.3", [], [], { "H2S|硫化氫|氫硫酸": {pg: "C2v", mp: "-85.5", bp: "-60.3", desc: "<strong>硫化氫</strong><br>具有腐敗雞蛋味的氣體，V型結構。", atoms: [{elem:"S",x:0,y:5,z:0,lpCount:2}, {elem:"H",x:40,y:-35,z:0}, {elem:"H",x:-40,y:-35,z:0}], bonds: [[0,1],[0,2]] }, "HS-|硫氫根": {pg: "Cinfv", mp: "-", bp: "-", desc: "<strong>氫硫根</strong><br>硫化氫的一級解離產物，硫原子有三對孤對電子。", atoms: [{elem:"S",x:-20,y:0,z:0,lpCount:3},{elem:"H",x:30,y:0,z:0}], bonds: [[0,1]] }});

// --- 4. 鹵化物系列 (全資料補完與鍵長修正) ---
const halideProps = { "BF3": ["-126.8", "-100.3"], "BCl3": ["-107", "12.6"], "BBr3": ["-46", "91.3"], "BI3": ["49.9", "210"], "AlF3": ["1290 (昇華)", "-"], "AlCl3": ["192.4", "120 (昇華)"], "AlBr3": ["97.5", "255"], "AlI3": ["191", "360"], "CF4": ["-183.6", "-127.8"], "CCl4": ["-22.9", "76.7"], "CBr4": ["90.1", "189.5"], "CI4": ["171 (分解)", "-"], "SiF4": ["-90", "-86 (昇華)"], "SiCl4": ["-70", "57.7"], "SiBr4": ["5", "154"], "SiI4": ["120.5", "287.5"], "NF3": ["-206.8", "-129"], "NCl3": ["-40", "71"], "NBr3": ["-100", "爆炸"], "NI3": ["-", "爆炸"], "PF3": ["-151.5", "-101.8"], "PCl3": ["-93.6", "76.1"], "PBr3": ["-41.5", "173.2"], "PI3": ["61", "分解"], "OF2": ["-223.8", "-144.8"], "OCl2": ["-135", "2.0"], "OBr2": ["-", "-"], "OI2": ["-", "-"], "SF2": ["-", "-"], "SCl2": ["-121", "59 (分解)"], "SBr2": ["-", "-"], "SI2": ["-", "-"] };
const haloNames = {'F':'氟', 'Cl':'氯', 'Br':'溴', 'I':'碘'};

['F','Cl','Br','I'].forEach(X => {
    const hn = haloNames[X]; let rX = (X==='F'?35: (X==='Cl'?40: (X==='Br'?45:50)));
    
    // BX3 系列
    let p = halideProps[`B${X}3`] || ["-","-"]; 
    addMol(`B${X}3|三${hn}化硼`, "B", "sp²", ["平面三角形","Trigonal Planar"], "120°", p[0], p[1], getTrigPlanar("B", X, 35+rX), [[0,1],[0,2],[0,3]], null, null, "D3h");
    
    // AlX3 系列
    p = halideProps[`Al${X}3`] || ["-","-"]; 
    addMol(`Al${X}3|三${hn}化鋁`, "Al", "sp²", ["平面三角形","Trigonal Planar"], "120°", p[0], p[1], getTrigPlanar("Al", X, 40+rX), [[0,1],[0,2],[0,3]], null, null, "D3h");
    
    // CX4 系列
    p = halideProps[`C${X}4`] || ["-","-"]; 
    addMol(`C${X}4|四${hn}化碳|四${hn}甲烷`, "C", "sp³", ["四面體","Tetrahedral"], "109.5°", p[0], p[1], getTetra("C", X, 35+rX), [[0,1],[0,2],[0,3],[0,4]], null, null, "Td");
    
    // SiX4 系列
    if(X !== 'Cl') { 
        p = halideProps[`Si${X}4`] || ["-","-"]; 
        addMol(`Si${X}4|四${hn}化矽`, "Si", "sp³", ["四面體","Tetrahedral"], "109.5°", p[0], p[1], getTetra("Si", X, 40+rX), [[0,1],[0,2],[0,3],[0,4]], null, null, "Td"); 
    }
    
    // NX3 系列
    p = halideProps[`N${X}3`] || ["-","-"]; 
    let dN = 35+rX, hN=dN*0.85, vN=dN*0.5; 
    addMol(`N${X}3|三${hn}化氮`, "N", "sp³", ["角錐形","Pyramidal"], (X==='F'?"102.3°":(X==='Cl'?"107.1°":(X==='Br'?"108°":"110°"))), p[0], p[1], [{elem:"N",x:0,y:15,z:0,lp3d:[{x:0,y:1,z:0}]},{elem:X,x:0,y:-10,z:hN},{elem:X,x:hN*0.866,y:-10,z:-hN*0.5},{elem:X,x:-hN*0.866,y:-10,z:-hN*0.5}], [[0,1],[0,2],[0,3]], null, null, "C3v");
    
    // PX3 系列
    if(X !== 'Cl') { 
        p = halideProps[`P${X}3`] || ["-","-"]; 
        let dP = 40+rX, hP=dP*0.85, vP=dP*0.5; 
        addMol(`P${X}3|三${hn}化磷`, "P", "sp³", ["角錐形","Pyramidal"], (X==='F'?"97.8°":(X==='Cl'?"100.3°":(X==='Br'?"101.5°":"102°"))), p[0], p[1], [{elem:"P",x:0,y:20,z:0,lp3d:[{x:0,y:1,z:0}]},{elem:X,x:0,y:-15,z:hP},{elem:X,x:hP*0.866,y:-15,z:-hP*0.5},{elem:X,x:-hP*0.866,y:-15,z:-hP*0.5}], [[0,1],[0,2],[0,3]], null, null, "C3v"); 
    }
    
    // OX2 系列
    p = halideProps[`O${X}2`] || ["-","-"]; 
    let dO = 35+rX; 
    addMol(`O${X}2|二${hn}化氧`, "O", "sp³", ["角形","Bent"], (X==='F'?"103.3°":(X==='Cl'?"110.9°":"114°")), p[0], p[1], [{elem:"O",x:0,y:0,z:0,lpCount:2},{elem:X,x:dO*0.8,y:-dO*0.6,z:0},{elem:X,x:-dO*0.8,y:-dO*0.6,z:0}], [[0,1],[0,2]], null, null, "C2v");
    
    // SX2 系列
    p = halideProps[`S${X}2`] || ["-","-"]; 
    let dS = 40+rX; 
    addMol(`S${X}2|二${hn}化硫`, "S", "sp³", ["角形","Bent"], (X==='F'?"98.2°":(X==='Cl'?"102.7°":"104°")), p[0], p[1], [{elem:"S",x:0,y:0,z:0,lpCount:2},{elem:X,x:dS*0.85,y:-dS*0.55,z:0},{elem:X,x:-dS*0.85,y:-dS*0.55,z:0}], [[0,1],[0,2]], null, null, "C2v");
});


// [保留] SiCl4 詳細資料
addMol("SiCl4|四氯化矽|Silicon Tetrachloride", "Si", "sp³", ["四面體","Tetrahedral"], "109.5°", "-70", "57.7", getTetra("Si", "Cl", 80), [[0,1],[0,2],[0,3],[0,4]], null,
    `<div class="info-section">
        <div class="info-title">⚗️ 物質性質</div>
        <div class="info-body">
            <span class="highlight-title">1. 立體結構：</span>中心矽原子採取 <strong>sp³ 混成</strong>，與同族的四氯化碳 (CCl₄) 具有相同的<strong>正四面體</strong>幾何結構，鍵角為 <strong>109.5°</strong>。<br>
            <span class="highlight-title">2. 物理性質：</span>常溫下為無色、易揮發的液體，具有強烈的刺鼻氣味。雖然 Si-Cl 鍵是極性共價鍵，但由於分子對稱性高，偶極矩互相抵銷，整體為<strong>非極性分子</strong>。<br>
            <span class="highlight-title">3. 化學性質：</span>與化學性質安定的 CCl₄ 不同，SiCl₄ 極易發生<strong>水解反應</strong>。這是因為矽原子的原子半徑較大，且擁有<strong>空 d 軌域</strong>，能接受水分子的氧原子進行親核攻擊，反應後生成矽酸並產生大量的氯化氫 (HCl) 白煙。
        </div>
    </div>
    <div class="info-section" style="margin-top: 12px; border-top: 1px dashed rgba(255,255,255,0.2); padding-top: 10px;">
        <div class="info-title">🏭 生活應用</div>
        <div class="info-body">
            <span class="highlight-title">1. 晶片製造 (多晶矽)：</span>它是半導體產業的基石。透過<strong>西門子法 (Siemens process)</strong>，將高純度的 SiCl₄ 與氫氣在 1100°C 高溫下反應還原，可製造出純度高達 99.9999999% (9N) 的<strong>電子級多晶矽</strong>，用於生產電腦晶片與太陽能電池。<br>
            <span class="highlight-title">2. 光纖通訊核心：</span>在光纖製程中，SiCl₄ 是最關鍵的原料。透過氣相沉積法將其高溫氧化，能生成折射率極高且無雜質的二氧化矽 (SiO₂)，構成光纖內層傳輸訊號的玻璃核心。<br>
            <span class="highlight-title">3. 軍事煙霧彈：</span>早期軍事上利用其「極易水解」的特性製作煙霧彈。當液態 SiCl₄ 炸開接觸空氣中的水氣時，會瞬間產生極濃密的白色酸霧 (HCl)，能有效遮蔽視線，但因具有毒性與腐蝕性，現代已較少使用。
        </div>
    </div>`
,"Td");

// [保留] PCl3 詳細資料
addMol("PCl3|三氯化磷|Phosphorus Trichloride", "P", "sp³", ["角錐形","Pyramidal"], "96-100°", "-93.6", "76.1", [{elem:"P",x:0,y:20,z:0,lp3d:[{x:0,y:1,z:0}]},{elem:"Cl",x:0,y:-15,z:68},{elem:"Cl",x:59,y:-15,z:-34},{elem:"Cl",x:-59,y:-15,z:-34}], [[0,1],[0,2],[0,3]], null,
    `<div class="info-section">
        <div class="info-title">⚗️ 物質性質</div>
        <div class="info-body">
            <span class="highlight-title">1. 立體結構：</span>中心磷原子採取 <strong>sp³ 混成</strong>。由於具有一對未共用電子對 (Lone Pair)，其對鍵結電子的斥力較大，導致 P-Cl 鍵角被壓縮至約 <strong>100°</strong>，形成<strong>三角錐形</strong>結構。<br>
            <span class="highlight-title">2. 物理性質：</span>常溫下為無色或微黃色的液體，會發煙。具有較低的沸點與強烈刺鼻味，可溶於苯、氯仿等有機溶劑。<br>
            <span class="highlight-title">3. 化學性質：</span>P-Cl 鍵極性大且反應性極高，遇水會劇烈<strong>水解</strong>並放熱，生成亞磷酸 (H₃PO₃) 與鹽酸霧。因磷原子上有一對孤對電子，可作為<strong>路易斯鹼</strong>參與配位反應。
        </div>
    </div>
    <div class="info-section" style="margin-top: 12px; border-top: 1px dashed rgba(255,255,255,0.2); padding-top: 10px;">
        <div class="info-title">🏭 生活應用</div>
        <div class="info-body">
            <span class="highlight-title">1. 除草劑原料 (嘉磷塞)：</span>工業上最大宗的用途是作為中間體，用於合成廣效性除草劑<strong>嘉磷塞 (Glyphosate)</strong>，這是目前全球農業使用量最大的農藥之一。<br>
            <span class="highlight-title">2. 有機合成 (氯化劑)：</span>在製藥與有機化學實驗室中，它是不可或缺的試劑。專門用來將有機分子中的<strong>羥基 (-OH)</strong> 取代為氯原子，或是將羧酸轉化為活性極高的醯氯，是合成染料與藥物的重要步驟。<br>
            <span class="highlight-title">3. 塑膠添加劑：</span>可用於製造含磷的<strong>阻燃劑</strong>與塑化劑。這些添加劑能讓電子產品的塑膠外殼在受熱時不易燃燒，大幅提升產品安全性。
        </div>
    </div>`
,"C3v");




// --- 5. 碳與其他氧化物 (直線型/平面型) ---
addMol("CO2|二氧化碳|乾冰", "C", "sp", ["直線型","Linear"], "180°", "-78.5 (昇華)", "-56.6", getLinear("C","O", 70), [[0,1,"double"],[0,2,"double"]], null, null, "Dinfh");
addMol("CS2|二硫化碳", "C", "sp", ["直線型","Linear"], "180°", "-111.6", "46.2", getLinear("C","S", 75), [[0,1,"double"],[0,2,"double"]], null, null, "Dinfh");
addMol("BeCl2|二氯化鈹", "Be", "sp", ["直線型","Linear"], "180°", "399", "482", getLinear("Be","Cl", 75), [[0,1], [0,2]], null, null, "Dinfh");
addMol("BCl3|三氯化硼", "B", "sp²", ["平面三角形","Trigonal Planar"], "120°", "-107", "12.6", getTrigPlanar("B","Cl", 75), [[0,1], [0,2], [0,3]], null, null, "D3h");
addMol("SO2|二氧化硫", "S", "sp²", ["角形","Bent"], "119°", "-72", "-10", 
    [
        {elem:"S", x:0, y:15, z:0, lpCount:1, lp3d:[{x:0,y:1,z:0}]}, 
        {elem:"O", x:55, y:-30, z:0}, 
        {elem:"O", x:-55, y:-30, z:0}
    ], 
    // 預設給兩個雙鍵 (擴大八隅體狀態)，讓程式去切換
    [[0,1,"double"], [0,2,"double"]], null, null, "C2v");
addMol("SO3|三氧化硫", "S", "sp²", ["平面三角形","Trigonal Planar"], "120°", "16.9", "44.8", getTrigPlanar("S","O", 68), [[0,1,"double"],[0,2,"double"],[0,3,"double"]], null, null, "D3h");
addMol("O3|臭氧", "O", "sp²", ["角形","Bent"], "117°", "-192.2", "-112", [{elem:"O",x:0,y:10,z:0,lp3d:[{x:0,y:1,z:0}]},{elem:"O",x:60,y:-35,z:0},{elem:"O",x:-60,y:-35,z:0,lpCount:3}], [[0,1,"double"],[0,2,"coordinate"]], null, null, "C2v");
addMol("NO2|二氧化氮", "N", "sp²", ["角形","Bent"], "134°", "-11.2", "21.2", [{elem:"N",x:0,y:10,z:0,lp3d:[{x:0,y:1,z:0}],radical:true},{elem:"O",x:60,y:-35,z:0},{elem:"O",x:-60,y:-35,z:0,lpCount:3}], [[0,1,"double"],[0,2,"coordinate"]], null, null, "C2v");
addMol("N2O|一氧化二氮|笑氣", "N", "sp", ["直線型","Linear"], "180°", "-90.8", "-88.5", [{elem:"N",x:0,y:0,z:0,lpCount:0},{elem:"N",x:-65,y:0,z:0,lp3d:[{x:-1,y:0,z:0}]},{elem:"O",x:65,y:0,z:0}], [[0,1,"triple"],[0,2,"coordinate"]], null, null, "Cinfv");
addMol("NO|一氧化氮", "雙原子", "sp²", ["直線型","Linear"], "-", "-164", "-152", [{elem:"N",x:-32,y:0,z:0,radical:true,lp3d:[{x:-1,y:1,z:0},{x:-1,y:-1,z:0.4},{x:-1,y:-1,z:-0.4}]},{elem:"O",x:32,y:0,z:0}], [[0,1,"double"]], null, null, "Cinfv");

// --- 6. 離子與特殊無機分子 (含共振結構) ---
// --- SCN- 共振結構展示 (修改：以 N=C=S 為預設) ---
addMol("SCN-|硫氰酸根", "C", "sp", ["直線型","Linear"], "180°", "-", "-", [], [], {
    "SCN-|主要共振結構 (N=C=S)": {pg: "Cinfv", mp: "-", bp: "-", atoms: [{elem:"C", x:0, y:0, z:0}, {elem:"N", x:-65, y:0, z:0, lpCount:2},{elem:"S", x:85, y:0, z:0, lpCount:2}], bonds: [[0,1,"double"], [0,2,"double"]] },
    "SCN-|次要共振結構 (N≡C-S)": {pg: "Cinfv", mp: "-", bp: "-", atoms: [{elem:"C", x:0, y:0, z:0}, {elem:"N", x:-60, y:0, z:0, lpCount:1}, {elem:"S", x:90, y:0, z:0, lpCount:3}], bonds: [[0,1,"triple"], [0,2,"single"]] }
});
addMol("NO+|亞硝鎓離子", "N", "sp", ["直線型","Linear"], "180°", "-", "-", [{elem:"N",x:-30,y:0,z:0,lpCount:1}, {elem:"O",x:30,y:0,z:0,lpCount:1}], [[0,1,"triple"]], null, null, "Cinfv");
addMol("NO2+|硝鎓離子", "N", "sp", ["直線型","Linear"], "180°", "-", "-", [{elem:"N",x:0,y:0,z:0}, {elem:"O",x:-65,y:0,z:0}, {elem:"O",x:65,y:0,z:0}], [[0,1,"double"],[0,2,"double"]], null, null, "Dinfh");
addMol("N3-|疊氮酸根", "N", "sp", ["直線型","Linear"], "180°", "-", "-", [], [], {
    "N3-|主要共振結構 (N=N=N)": {pg: "Dinfh", atoms: [{elem:"N",x:0,y:0,z:0},{elem:"N",x:-65,y:0,z:0,lpCount:2},{elem:"N",x:65,y:0,z:0,lpCount:2}], bonds: [[0,1,"double"],[0,2,"double"]] },
    "N3-|主要共振結構 (N≡N-N)": {pg: "Dinfh", atoms: [{elem:"N",x:0,y:0,z:0},{elem:"N",x:-60,y:0,z:0,lpCount:1},{elem:"N",x:85,y:0,z:0,lpCount:3}], bonds: [[0,1,"triple"],[0,2,"single"]] }
});
addMol("OCN-|氰酸根", "C", "sp", ["直線型","Linear"], "180°", "-", "-", [], [], {
    "OCN-|主要共振結構 (N≡C-O)": {pg: "Cinfv", atoms: [{elem:"C",x:0,y:0,z:0},{elem:"N",x:-60,y:0,z:0,lpCount:1},{elem:"O",x:85,y:0,z:0,lpCount:3}], bonds: [[0,1,"triple"],[0,2,"single"]] },
    "OCN-|次要共振結構 (N=C=O)": {pg: "Cinfv", atoms: [{elem:"C",x:0,y:0,z:0},{elem:"N",x:-65,y:0,z:0,lpCount:2},{elem:"O",x:65,y:0,z:0,lpCount:2}], bonds: [[0,1,"double"],[0,2,"double"]] }
});
addMol("CNO-|雷酸根", "N", "sp", ["直線型","Linear"], "180°", "-", "-", [], [], {
    "CNO-|主要共振結構 (C≡N-O)": {pg: "Cinfv", atoms: [{elem:"N",x:0,y:0,z:0},{elem:"C",x:-60,y:0,z:0,lpCount:1},{elem:"O",x:85,y:0,z:0,lpCount:3}], bonds: [[0,1,"triple"],[0,2,"single"]] },
    "CNO-|次要共振結構 (C=N=O)": {pg: "Cinfv", atoms: [{elem:"N",x:0,y:0,z:0},{elem:"C",x:-65,y:0,z:0,lpCount:2},{elem:"O",x:65,y:0,z:0,lpCount:2}], bonds: [[0,1,"double"],[0,2,"double"]] }
});
addMol("HOCN|氰酸", "C", "sp", ["直線/角形","Linear/Bent"], "180°/105°", "-86", "23.5", [{elem:"C",x:0,y:0,z:0}, {elem:"N",x:65,y:0,z:0,lpCount:1}, {elem:"O",x:-65,y:0,z:0,lpCount:2}, {elem:"H",x:-95,y:30,z:0}], [[0,1,"triple"], [0,2], [2,3]], null, null, "Cs");

// --- 7. 擴大八隅體與複雜幾何構型 ---
// P-Cl=80, P-Br=85, S-F=75, S=O=68, Xe=O=76
addMol("PCl5|五氯化磷", "P", "sp³d", ["雙三角錐","Trigonal Bipyramidal"], "90°, 120°", "160.5", "166.8", [{elem:"P",x:0,y:0,z:0},{elem:"Cl",x:0,y:0,z:85},{elem:"Cl",x:0,y:0,z:-85},{elem:"Cl",x:-80,y:0,z:0},{elem:"Cl",x:40,y:69,z:0},{elem:"Cl",x:40,y:-69,z:0}], [[0,1],[0,2],[0,3],[0,4],[0,5]], null, null, "D3h");
addMol("PBr5|五溴化磷", "P", "sp³d", ["雙三角錐","Trigonal Bipyramidal"], "90°, 120°", "100 (分解)", "106 (分解)", [{elem:"P",x:0,y:0,z:0}, {elem:"Br",x:0,y:90,z:0}, {elem:"Br",x:0,y:-90,z:0}, {elem:"Br",x:85,y:0,z:0}, {elem:"Br",x:-42,y:0,z:74}, {elem:"Br",x:-42,y:0,z:-74}], [[0,1],[0,2],[0,3],[0,4],[0,5]], null, null, "D3h");
addMol("SF6|六氟化硫", "S", "sp³d²", ["八面體","Octahedral"], "90°", "-50.8", "-63.8 (昇華)", getOcta("S","F", 75), [[0,1],[0,2],[0,3],[0,4],[0,5],[0,6]], null, null, "Oh");
addMol("SF4|四氟化硫", "S", "sp³d", ["翹翹板型","Seesaw"], "<90°, <120°", "-121", "-38", [{elem:"S",x:0,y:0,z:0,lp3d:[{x:0,y:-1,z:0}]},{elem:"F",x:0,y:0,z:80},{elem:"F",x:0,y:0,z:-80},{elem:"F",x:45,y:65,z:0},{elem:"F",x:-45,y:65,z:0}], [[0,1],[0,2],[0,3],[0,4]], null, null, "C2v");
addMol("ClF3|三氟化氯", "Cl", "sp³d", ["T型","T-shaped"], "<90°", "-76.3", "11.8", [{elem:"Cl",x:0,y:0,z:0,lp3d:[{x:-1,y:0.5,z:0}, {x:-1,y:-0.5,z:0}]}, {elem:"F",x:0,y:80,z:0}, {elem:"F",x:0,y:-80,z:0}, {elem:"F",x:70,y:0,z:0}], [[0,1],[0,2],[0,3]], null, null, "C2v");
addMol("XeF2|二氟化氙", "Xe", "sp³d", ["直線型","Linear"], "180°", "128.6", "-", [{elem:"Xe",x:0,y:0,z:0,lp3d:[{x:0,y:1,z:0}, {x:0.866,y:-0.5,z:0}, {x:-0.866,y:-0.5,z:0}]}, {elem:"F",x:0,y:0,z:85}, {elem:"F",x:0,y:0,z:-85}], [[0,1],[0,2]], null, null, "Dinfh");

addMol("XeF4|四氟化氙", "Xe", "sp³d²", ["平面四邊形","Square Planar"], "90°", "117 (昇華)", "-", [
{elem:"Xe",x:0,y:0,z:0,lp3d:[{x:1,y:0,z:0}, {x:-1,y:0,z:0}]}, 
{elem:"F",x:0,y:85,z:0}, {elem:"F",x:0,y:-85,z:0}, 
{elem:"F",x:0,y:0,z:85}, {elem:"F",x:0,y:0,z:-85}
], [[0,1],[0,2],[0,3],[0,4]], null, null, "D4h");

addMol("BrF5|五氟化溴", "Br", "sp³d²", ["四角錐","Square Pyramidal"], "<90°", "-61.3", "40.3", [{elem:"Br",x:0,y:0,z:0,lp3d:[{x:0,y:-1,z:0}]}, {elem:"F",x:0,y:80,z:0}, {elem:"F",x:70,y:0,z:0}, {elem:"F",x:-70,y:0,z:0}, {elem:"F",x:0,y:0,z:70}, {elem:"F",x:0,y:0,z:-70}], [[0,1],[0,2],[0,3],[0,4],[0,5]], null, null, "C4v");
addMol("IF7|七氟化碘", "I", "sp³d³", ["五角雙錐","Pentagonal Bipyramidal"], "72°, 90°", "4.8", "4.8 (昇華)", [{elem:"I",x:0,y:0,z:0,lpCount:0}, {elem:"F",x:0,y:90,z:0}, {elem:"F",x:0,y:-90,z:0}, {elem:"F",x:80,y:0,z:0}, {elem:"F",x:25,y:0,z:76}, {elem:"F",x:25,y:0,z:-76}, {elem:"F",x:-65,y:0,z:47}, {elem:"F",x:-65,y:0,z:-47}], [[0,1],[0,2],[0,3],[0,4],[0,5],[0,6],[0,7]], null, null, "D5h");
addMol("SeF6|六氟化硒", "Se", "sp³d²", ["八面體","Octahedral"], "90°", "-34.6", "-46.6 (昇華)", getOcta("Se","F", 75), [[0,1],[0,2],[0,3],[0,4],[0,5],[0,6]], null, null, "Oh");
addMol("TeF6|六氟化碲", "Te", "sp³d²", ["八面體","Octahedral"], "90°", "-37.6", "-38.9 (昇華)", getOcta("Te","F", 75), [[0,1],[0,2],[0,3],[0,4],[0,5],[0,6]], null, null, "Oh");
addMol("AsF5|五氟化砷", "As", "sp³d", ["雙三角錐","Trigonal Bipyramidal"], "90°, 120°", "-79.8", "-52.8", [{elem:"As",x:0,y:0,z:0}, {elem:"F",x:0,y:80,z:0}, {elem:"F",x:0,y:-80,z:0}, {elem:"F",x:70,y:0,z:0}, {elem:"F",x:-35,y:0,z:60}, {elem:"F",x:-35,y:0,z:-60}], [[0,1],[0,2],[0,3],[0,4],[0,5]], null, null, "D3h");
addMol("TeF4|四氟化碲", "Te", "sp³d", ["翹翹板型","Seesaw"], "<90°, <120°", "129.6", "193", [{elem:"Te",x:0,y:0,z:0,lp3d:[{x:-1,y:0,z:0}]}, {elem:"F",x:0,y:85,z:0}, {elem:"F",x:0,y:-85,z:0}, {elem:"F",x:70,y:0,z:50}, {elem:"F",x:70,y:0,z:-50}], [[0,1],[0,2],[0,3],[0,4]], null, null, "C2v");
addMol("XeO3|三氧化氙", "Xe", "sp³", ["角錐形","Trigonal Pyramidal"], "103°", "25 (爆炸)", "-", [{elem:"Xe",x:0,y:20,z:0,lp3d:[{x:0,y:1,z:0}]},{elem:"O",x:0,y:-30,z:57},{elem:"O",x:49,y:-30,z:-28.5},{elem:"O",x:-49,y:-30,z:-28.5}], [[0,1,"double"],[0,2,"double"],[0,3,"double"]], null, null, "C3v");
addMol("XeO4|四氧化氙", "Xe", "sp³", ["四面體","Tetrahedral"], "109.5°", "-35.9", "0 (分解)", getTetra("Xe","O", 76), [[0,1,"double"],[0,2,"double"],[0,3,"double"],[0,4,"double"]], null, null, "Td");
addMol("XeOF4|四氟氧化氙|XeOF4", "Xe", "sp³d²", ["四角錐","Square Pyramidal"], "<90°", "-46", "101", [{elem:"Xe",x:0,y:0,z:0,lp3d:[{x:0,y:-1,z:0}]}, {elem:"O",x:0,y:80,z:0}, {elem:"F",x:80,y:0,z:0}, {elem:"F",x:-80,y:0,z:0}, {elem:"F",x:0,y:0,z:80}, {elem:"F",x:0,y:0,z:-80}], [[0,1,"double"],[0,2],[0,3],[0,4],[0,5]], null, null, "C4v");
addMol("IOF5|五氟氧化碘", "I", "sp³d²", ["八面體","Octahedral"], "90°", "4.5", "110", [{elem:"I",x:0,y:0,z:0}, {elem:"O",x:0,y:85,z:0}, {elem:"F",x:0,y:-85,z:0}, {elem:"F",x:85,y:0,z:0}, {elem:"F",x:-85,y:0,z:0}, {elem:"F",x:0,y:0,z:85}, {elem:"F",x:0,y:0,z:-85}], [[0,1,"double"],[0,2],[0,3],[0,4],[0,5],[0,6]], null, null, "C4v");
addMol("AsF3|三氟化砷", "As", "sp³", ["角錐形","Pyramidal"], "96°", "-6", "57.8", [{elem:"As",x:0,y:15,z:0,lp3d:[{x:0,y:1,z:0}]}, {elem:"F",x:0,y:-45,z:55}, {elem:"F",x:48,y:-45,z:-28}, {elem:"F",x:-48,y:-45,z:-28}], [[0,1],[0,2],[0,3]], null, null, "C3v");
addMol("SbCl3|三氯化銻", "Sb", "sp³", ["角錐形","Pyramidal"], "97°", "73.4", "220.3", [{elem:"Sb",x:0,y:15,z:0,lp3d:[{x:0,y:1,z:0}]}, {elem:"Cl",x:0,y:-55,z:65}, {elem:"Cl",x:55,y:-55,z:-35}, {elem:"Cl",x:-55,y:-55,z:-35}], [[0,1],[0,2],[0,3]], null, null, "C3v");
addMol("ICl3|三氯化碘|Iodine Trichloride", "I", "sp³d", ["T型","T-shaped"], "<90°", "101 (分解)", "-", [{elem:"I",x:0,y:0,z:0,lpCount:2,lp3d:[{x:-1,y:0.5,z:0},{x:-1,y:-0.5,z:0}]},{elem:"Cl",x:90,y:0,z:0},{elem:"Cl",x:0,y:90,z:0},{elem:"Cl",x:0,y:-90,z:0}], [[0,1],[0,2],[0,3]], null, '<div class="info-section"><div class="info-title">🧪 物質簡介</div><div class="info-body"><strong>三氯化碘 (ICl₃)</strong><br>中心碘原子採取 sp³d 混成。為了減少電子雲斥力，兩對孤對電子佔據水平位置，使分子呈現 T 型結構。</div></div>', "C2v");
addMol("B2H6|乙硼烷|Diborane", "B", "sp³", ["特殊 (含氫橋鍵)","Banana Bonds"], "120°(端)/97°(橋)", "-164.8", "-92.5", [{elem:"B",x:-40,y:0,z:0,lpCount:0},{elem:"B",x:40,y:0,z:0,lpCount:0},{elem:"H",x:0,y:0,z:50},{elem:"H",x:0,y:0,z:-50},{elem:"H",x:-65,y:43,z:0},{elem:"H",x:-65,y:-43,z:0},{elem:"H",x:65,y:43,z:0},{elem:"H",x:65,y:-43,z:0}], [[0,2],[0,3],[1,2],[1,3],[0,4],[0,5],[1,6],[1,7]], null, '<div class="info-section"><div class="info-title">🍌 結構特性</div><div class="info-body"><strong>乙硼烷 (B₂H₆)</strong><br>具有三中心二電子鍵。每個硼原子與四個氫原子連線，形成類似 sp³ 的幾何排列。</div></div>', "D2h");
// --- 8. 陰離子群 (Complex Anions) ---
addMol("SiF62-|六氟矽酸根", "Si", "sp³d²", ["八面體","Octahedral"], "90°", "-", "-", getOcta("Si","F", 75), [[0,1],[0,2],[0,3],[0,4],[0,5],[0,6]], null, null, "Oh");
addMol("PF6-|六氟磷酸根", "P", "sp³d²", ["八面體","Octahedral"], "90°", "-", "-", getOcta("P","F", 75), [[0,1],[0,2],[0,3],[0,4],[0,5],[0,6]], null, null, "Oh");
addMol("SbF6-|六氟銻酸根", "Sb", "sp³d²", ["八面體","Octahedral"], "90°", "-", "-", getOcta("Sb","F", 80), [[0,1],[0,2],[0,3],[0,4],[0,5],[0,6]], null, null, "Oh");
addMol("I3-|三碘陰離子|三碘錯離子", "I", "sp³d", ["直線型","Linear"], "180°", "-", "-", [{elem:"I",x:0,y:0,z:0,lp3d:[{x:0,y:1,z:0},{x:0,y:-0.5,z:0.866},{x:0,y:-0.5,z:-0.866}]},{elem:"I",x:-100,y:0,z:0},{elem:"I",x:100,y:0,z:0}], [[0,1],[0,2]], null, null, "Dinfh");
addMol("ICl2-|二氯碘離子", "I", "sp³d", ["直線型","Linear"], "180°", "-", "-", [{elem:"I",x:0,y:0,z:0,lp3d:[{x:0,y:1,z:0},{x:0,y:-0.5,z:0.866},{x:0,y:-0.5,z:-0.866}]}, {elem:"Cl",x:-90,y:0,z:0}, {elem:"Cl",x:90,y:0,z:0}], [[0,1],[0,2]], null, null, "Dinfh");
addMol("ICl4-|四氯碘離子", "I", "sp³d²", ["平面四邊形","Square Planar"], "90°", "-", "-", [{elem:"I",x:0,y:0,z:0,lp3d:[{x:0,y:1,z:0}, {x:0,y:-1,z:0}]}, {elem:"Cl",x:90,y:0,z:0}, {elem:"Cl",x:-90,y:0,z:0}, {elem:"Cl",x:0,y:0,z:90}, {elem:"Cl",x:0,y:0,z:-90}], [[0,1],[0,2],[0,3],[0,4]], null, null, "D4h");
addMol("BF4-|四氟硼酸根", "B", "sp³", ["四面體","Tetrahedral"], "109.5°", "-", "-", getTetra("B","F", 70), [[0,1],[0,2],[0,3],[0,4]], null, null, "Td");
addMol("BBF4-|四氟硼酸根", "B", "sp³", ["四面體","Tetrahedral"], "109.5°", "-", "-", getTetra("B","F", 70), [[0,1],[0,2],[0,3],[4,0,"coordinate"]], null, null, "Td");
addMol("AlCl4-|四氯鋁酸根", "Al", "sp³", ["四面體","Tetrahedral"], "109.5°", "-", "-", getTetra("Al","Cl", 80), [[0,1],[0,2],[0,3],[0,4]], null, null, "Td");
addMol("BH4-|硼氫化離子", "B", "sp³", ["四面體","Tetrahedral"], "109.5°", "-", "-", getTetra("B","H", 50), [[0,1],[0,2],[0,3],[0,4]], null, null, "Td");

// --- 9. 酸根與含氧酸 ---
// --- 酸根與含氧酸 (修正離子鍵距離與鍵級顯示) ---
addMol("H2SO4|硫酸系列", "S", "sp³", ["四面體","Tetrahedral"], "109.5°", "10.3", "337", [], [], {
    "H2SO4|硫酸": { pg: "C2", mp: "10.3", bp: "337", desc: "<strong>硫酸</strong><br>工業之母，具強脫水性與氧化性，由兩個配位鍵 (S→O) 與兩個 S-OH 構成，分子電中性。", atoms: [{elem:"S",x:0,y:0,z:0},{elem:"O",x:0,y:68,z:0,lpCount:3},{elem:"O",x:0,y:-25,z:-63,lpCount:3},{elem:"O",x:60,y:-30,z:35,lpCount:2},{elem:"O",x:-60,y:-30,z:35,lpCount:2},{elem:"H",x:85,y:5,z:60},{elem:"H",x:-85,y:5,z:60}], bonds: [[0,1,"coordinate"],[0,2,"coordinate"],[0,3,"single"],[0,4,"single"],[3,5],[4,6]] },
    "HSO4-|硫酸氫根": { pg: "Cs", mp: "-", bp: "-", desc: "<strong>硫酸氫根</strong><br>酸式鹽陰離子，水溶液呈強酸性，S-O⁻ 端帶有負電荷。", atoms: [{elem:"S",x:0,y:0,z:0},{elem:"O",x:0,y:68,z:0,lpCount:3},{elem:"O",x:0,y:-25,z:-63,lpCount:3},{elem:"O",x:60,y:-30,z:35,lpCount:3},{elem:"O",x:-60,y:-30,z:35,lpCount:2},{elem:"H",x:-85,y:5,z:60}], bonds: [[0,1,"coordinate"],[0,2,"coordinate"],[0,3,"single"],[0,4,"single"],[4,5]] },
    "SO42-|硫酸根": { pg: "Td", mp: "-", bp: "-", desc: "<strong>硫酸根</strong><br>正四面體結構，化學性質穩定，兩個 S-O⁻ 端顯示粉紅電子。", atoms: [{elem:"S",x:0,y:0,z:0},{elem:"O",x:0,y:68,z:0,lpCount:3},{elem:"O",x:0,y:-25,z:-63,lpCount:3},{elem:"O",x:60,y:-30,z:35,lpCount:3},{elem:"O",x:-60,y:-30,z:35,lpCount:3}], bonds: [[0,1,"coordinate"],[0,2,"coordinate"],[0,3,"single"],[0,4,"single"]] },
    "NaHSO4|硫酸氫鈉": { pg: "Cs", mp: "58 (分解)", bp: "-", desc: "<strong>硫酸氫鈉</strong><br>溶於水呈強酸性，常用於清潔劑或降低 pH 值。", atoms: [{elem:"S",x:-20,y:0,z:0},{elem:"O",x:-20,y:68,z:0,lpCount:3},{elem:"O",x:-20,y:-25,z:-63,lpCount:3},{elem:"O",x:40,y:-30,z:35,lpCount:3},{elem:"O",x:-80,y:-30,z:35,lpCount:2},{elem:"H",x:-105,y:-5,z:60},{elem:"Na",x:100,y:40,z:0,r:15}], bonds: [[0,1,"coordinate"],[0,2,"coordinate"],[0,3,"single"],[0,4,"single"],[4,5]] },
    "KHSO4|硫酸氫鉀": { pg: "Cs", mp: "197", bp: "-", desc: "<strong>硫酸氫鉀</strong><br>易溶於水呈強酸性，加熱失水可製備焦硫酸鉀。", atoms: [{elem:"S",x:-20,y:0,z:0},{elem:"O",x:-20,y:68,z:0,lpCount:3},{elem:"O",x:-20,y:-25,z:-63,lpCount:3},{elem:"O",x:40,y:-30,z:35,lpCount:3},{elem:"O",x:-80,y:-30,z:35,lpCount:2},{elem:"H",x:-105,y:-5,z:60},{elem:"K",x:110,y:40,z:0,r:22}], bonds: [[0,1,"coordinate"],[0,2,"coordinate"],[0,3,"single"],[0,4,"single"],[4,5]] },
    "CaSO4|硫酸鈣|石膏": { pg: "Td", mp: "1460", bp: "-", desc: "<strong>硫酸鈣 (石膏)</strong><br>微溶於水，廣泛用於建築材料、模型製作與作為豆腐凝固劑。", atoms: [{elem:"S",x:0,y:0,z:0},{elem:"O",x:0,y:68,z:0,lpCount:3},{elem:"O",x:0,y:-25,z:-63,lpCount:3},{elem:"O",x:60,y:-30,z:35,lpCount:3},{elem:"O",x:-60,y:-30,z:35,lpCount:3},{elem:"Ca",x:0,y:0,z:100,r:20,lpCount:0}], bonds: [[0,1,"coordinate"],[0,2,"coordinate"],[0,3,"single"],[0,4,"single"]] },
    "BaSO4|硫酸鋇|重晶石": { pg: "Td", mp: "1580", bp: "-", desc: "<strong>硫酸鋇 (重晶石)</strong><br>極難溶於水與酸，無毒且密度大，醫學上用於消化道X光攝影(鋇餐)。", atoms: [{elem:"S",x:0,y:0,z:0},{elem:"O",x:0,y:68,z:0,lpCount:3},{elem:"O",x:0,y:-25,z:-63,lpCount:3},{elem:"O",x:60,y:-30,z:35,lpCount:3},{elem:"O",x:-60,y:-30,z:35,lpCount:3},{elem:"Ba",x:0,y:0,z:110,r:25,lpCount:0}], bonds: [[0,1,"coordinate"],[0,2,"coordinate"],[0,3,"single"],[0,4,"single"]] },
    "CuSO4|硫酸銅": { pg: "Td", mp: "110 (失水)", bp: "-", desc: "<strong>硫酸銅</strong><br>無水物為白色，吸水後變藍色(五水合)，常用於游泳池殺菌、波爾多液原料與電鍍。", atoms: [{elem:"S",x:0,y:0,z:0},{elem:"O",x:0,y:68,z:0,lpCount:3},{elem:"O",x:0,y:-25,z:-63,lpCount:3},{elem:"O",x:60,y:-30,z:35,lpCount:3},{elem:"O",x:-60,y:-30,z:35,lpCount:3},{elem:"Cu",x:0,y:0,z:100,r:18,lpCount:0}], bonds: [[0,1,"coordinate"],[0,2,"coordinate"],[0,3,"single"],[0,4,"single"]] },
    "FeSO4|硫酸亞鐵|綠礬": { pg: "Td", mp: "64 (失水)", bp: "-", desc: "<strong>硫酸亞鐵 (綠礬)</strong><br>淺綠色晶體，常用於醫療補血劑(鐵劑)、水處理絮凝劑與還原劑。", atoms: [{elem:"S",x:0,y:0,z:0},{elem:"O",x:0,y:68,z:0,lpCount:3},{elem:"O",x:0,y:-25,z:-63,lpCount:3},{elem:"O",x:60,y:-30,z:35,lpCount:3},{elem:"O",x:-60,y:-30,z:35,lpCount:3},{elem:"Fe",x:0,y:0,z:100,r:18,lpCount:0}], bonds: [[0,1,"coordinate"],[0,2,"coordinate"],[0,3,"single"],[0,4,"single"]] },
    "ZnSO4|硫酸鋅|皓礬": { pg: "Td", mp: "100 (失水)", bp: "500 (分解)", desc: "<strong>硫酸鋅 (皓礬)</strong><br>無色針狀晶體，用於製造人造纖維、木材防腐與農業微量元素肥料。", atoms: [{elem:"S",x:0,y:0,z:0},{elem:"O",x:0,y:68,z:0,lpCount:3},{elem:"O",x:0,y:-25,z:-63,lpCount:3},{elem:"O",x:60,y:-30,z:35,lpCount:3},{elem:"O",x:-60,y:-30,z:35,lpCount:3},{elem:"Zn",x:0,y:0,z:100,r:18,lpCount:0}], bonds: [[0,1,"coordinate"],[0,2,"coordinate"],[0,3,"single"],[0,4,"single"]] },
    "MgSO4|硫酸鎂|瀉鹽": { pg: "Td", mp: "1124", bp: "-", desc: "<strong>硫酸鎂 (瀉鹽)</strong><br>易溶於水，醫療上作為瀉劑或緩解子癇，生活中常用於泡澡浴鹽放鬆肌肉。", atoms: [{elem:"S",x:0,y:0,z:0},{elem:"O",x:0,y:68,z:0,lpCount:3},{elem:"O",x:0,y:-25,z:-63,lpCount:3},{elem:"O",x:60,y:-30,z:35,lpCount:3},{elem:"O",x:-60,y:-30,z:35,lpCount:3},{elem:"Mg",x:0,y:0,z:100,r:18,lpCount:0}], bonds: [[0,1,"coordinate"],[0,2,"coordinate"],[0,3,"single"],[0,4,"single"]] }
});

addMol("H2SO3|亞硫酸系列", "S", "sp³", ["角錐形","Pyramidal"], "106°", "-", "不穩定", [], [], {
"H2SO3|亞硫酸": { pg: "Cs", mp: "-", bp: "不穩定", desc: "<strong>亞硫酸</strong><br>僅存在於水溶液中的二元弱酸，極不穩定。具有強還原性與漂白能力，受熱或久置易分解出二氧化硫氣體。", atoms: [{elem:"S",x:0,y:15,z:0,lpCount:1},{elem:"O",x:0,y:80,z:0},{elem:"O",x:55,y:-30,z:30},{elem:"O",x:-55,y:-30,z:30},{elem:"H",x:85,y:-10,z:30},{elem:"H",x:-85,y:-10,z:30}], bonds: [[0,1,"double"],[0,2],[0,3],[2,4],[3,5]] },
"HSO3-|亞硫酸氫根": { pg: "Cs", mp: "-", bp: "-", desc: "<strong>亞硫酸氫根</strong><br>亞硫酸的第一級電離產物，為兩性離子。在酸性環境中不穩定，廣泛存在於亞硫酸氫鹽溶液中，具抗氧化性質。", atoms: [{elem:"S",x:0,y:15,z:0,lpCount:1},{elem:"O",x:0,y:80,z:0},{elem:"O",x:55,y:-30,z:30},{elem:"O",x:-55,y:-30,z:30},{elem:"H",x:85,y:-10,z:30}], bonds: [[0,1,"double"],[0,2],[0,3],[2,4]] },
"SO32-|亞硫酸根": { pg: "C3v", mp: "-", bp: "-", desc: "<strong>亞硫酸根</strong><br>亞硫酸的完全電離產物，中心硫原子有一對孤對電子。具有強還原性，易被空氣中的氧氧化成硫酸根。", atoms: [{elem:"S",x:0,y:15,z:0,lpCount:1},{elem:"O",x:0,y:80,z:0},{elem:"O",x:55,y:-30,z:30},{elem:"O",x:-55,y:-30,z:30}], bonds: [[0,1,"double"],[0,2,"single"],[0,3,"single"]] },
"Na2SO3|亞硫酸鈉": { pg: "C3v", mp: "33.4 (分解)", bp: "-", desc: "<strong>亞硫酸鈉</strong><br>常見的亞硫酸鹽，為白色粉末，易溶於水。常用作還原劑、防腐劑以及攝影顯影劑的保護劑。", atoms: [{elem:"S",x:0,y:15,z:0,lpCount:1},{elem:"O",x:0,y:80,z:0},{elem:"O",x:55,y:-30,z:30},{elem:"O",x:-55,y:-30,z:30},{elem:"Na",x:100,y:20,z:0,r:15},{elem:"Na",x:-100,y:20,z:0,r:15}], bonds: [[0,1,"double"],[0,2,"single"],[0,3,"single"],[4,2,"ionic_thin"],[5,3,"ionic_thin"]] },
"NaHSO3|亞硫酸氫鈉": { pg: "Cs", mp: "150 (分解)", bp: "-", desc: "<strong>亞硫酸氫鈉</strong><br>亞硫酸的酸式鹽，為白色結晶粉末，有二氧化硫的刺激氣氣味。常用於漂白織物、食品防腐及處理工業廢水。", atoms: [{elem:"S",x:0,y:15,z:0,lpCount:1},{elem:"O",x:0,y:80,z:0},{elem:"O",x:55,y:-30,z:30},{elem:"O",x:-55,y:-30,z:30},{elem:"H",x:85,y:-10,z:30},{elem:"Na",x:-100,y:0,z:0,r:15}], bonds: [[0,1,"double"],[0,2],[0,3],[2,4],[5,3,"ionic_thin"]] },
"CaSO3|亞硫酸鈣": { pg: "C3v", mp: "600 (分解)", bp: "-", desc: "<strong>亞硫酸鈣</strong><br>白色結晶粉末，微溶於水。主要用作食品防腐劑、消毒劑，也是煙氣脫硫工藝中的常見產物。", atoms: [{elem:"S",x:0,y:15,z:0,lpCount:1},{elem:"O",x:0,y:80,z:0},{elem:"O",x:55,y:-30,z:30},{elem:"O",x:-55,y:-30,z:30},{elem:"Ca",x:0,y:0,z:90,r:20}], bonds: [[0,1,"double"],[0,2,"single"],[0,3,"single"],[4,2,"ionic_thin"],[4,3,"ionic_thin"]] }
});

addMol("H2S2O3|硫代硫酸系列", "S", "sp³", ["四面體","Tetrahedral"], "109.5°", "-78 (分解)", "-", [], [], {
    "H2S2O3|硫代硫酸": { pg: "Cs", mp: "-78", bp: "-", desc: "<strong>硫代硫酸</strong><br>不穩定酸，中心S連接另一個外圍S原子。", atoms: [{elem:"S",x:0,y:0,z:0},{elem:"S",x:0,y:80,z:0},{elem:"O",x:0,y:-25,z:-63},{elem:"O",x:60,y:-30,z:35},{elem:"O",x:-60,y:-30,z:35},{elem:"H",x:85,y:5,z:60},{elem:"H",x:-85,y:5,z:60}], bonds: [[0,1,"double"],[0,2,"double"],[0,3],[0,4],[3,5],[4,6]] },
    "HS2O3-|硫代硫酸氫根": { pg: "Cs", mp: "-", bp: "-", desc: "<strong>硫代硫酸氫根</strong><br>結構類似硫酸氫根但一個O被S取代。", atoms: [{elem:"S",x:0,y:0,z:0},{elem:"S",x:0,y:80,z:0},{elem:"O",x:0,y:-25,z:-63},{elem:"O",x:60,y:-30,z:35},{elem:"O",x:-60,y:-30,z:35},{elem:"H",x:85,y:5,z:60}], bonds: [[0,1,"double"],[0,2,"double"],[0,3],[0,4],[3,5]] },
    "S2O32-|硫代硫酸根": { pg: "C3v", mp: "-", bp: "-", desc: "<strong>硫代硫酸根</strong><br>具還原性，中心硫原子與外圍硫形成雙鍵。", atoms: [{elem:"S",x:0,y:0,z:0},{elem:"S",x:0,y:80,z:0},{elem:"O",x:0,y:-25,z:-63},{elem:"O",x:60,y:-30,z:35},{elem:"O",x:-60,y:-30,z:35}], bonds: [[0,1,"double"],[0,2,"double"],[0,3,"single"],[0,4,"single"]] },
    "Na2S2O3|硫代硫酸鈉|大蘇打|海波": { pg: "C3v", mp: "48.3", bp: "100 (分解)", desc: "<strong>硫代硫酸鈉 (海波)</strong><br>Na⁺ 位於結構外側，無實體鍵連線。", atoms: [{elem:"S",x:0,y:0,z:0},{elem:"S",x:0,y:80,z:0},{elem:"O",x:0,y:-25,z:-63},{elem:"O",x:60,y:-30,z:35},{elem:"O",x:-60,y:-30,z:35},{elem:"Na",x:100,y:20,z:0,r:15},{elem:"Na",x:-100,y:20,z:0,r:15}], bonds: [[0,1,"double"],[0,2,"double"],[0,3,"single"],[0,4,"single"]] }
});

addMol("H2CO3|碳酸系列", "C", "sp²", ["平面三角形","Trigonal Planar"], "120°", "-", "不穩定", [], [], {
    "H2CO3|碳酸": { pg: "Cs", mp: "-", bp: "不穩定", desc: "<strong>碳酸</strong><br>二質子弱酸，存在於汽水中。", atoms: [{elem:"C",x:0,y:0,z:0},{elem:"O",x:0,y:70,z:0},{elem:"O",x:60,y:-35,z:0},{elem:"O",x:-60,y:-35,z:0},{elem:"H",x:90,y:-10,z:0},{elem:"H",x:-90,y:-10,z:0}], bonds: [[0,1,"double"],[0,2],[0,3],[2,4],[3,5]] },
    "HCO3-|碳酸氫根": { pg: "Cs", mp: "-", bp: "-", desc: "<strong>碳酸氫根</strong><br>帶-1價電荷，小蘇打的主要成分。", atoms: [{elem:"C",x:0,y:0,z:0},{elem:"O",x:0,y:70,z:0},{elem:"O",x:60,y:-35,z:0},{elem:"O",x:-60,y:-35,z:0},{elem:"H",x:-90,y:-10,z:0}], bonds: [[0,1,"double"],[0,2],[0,3],[3,4]] },
    "CO32-|碳酸根": { pg: "D3h", mp: "-", bp: "-", desc: "<strong>碳酸根</strong><br>帶-2價電荷，共振結構。", atoms: [{elem:"C",x:0,y:0,z:0},{elem:"O",x:0,y:70,z:0},{elem:"O",x:60,y:-35,z:0},{elem:"O",x:-60,y:-35,z:0}], bonds: [[0,1,"double"],[0,2,"single"],[0,3,"single"]] },
    "CaCO3|碳酸鈣|灰石": { pg: "D3h", mp: "825 (分解)", bp: "-", desc: "<strong>碳酸鈣</strong><br>Ca²⁺ 位於碳酸根平面上方。", atoms: [{elem:"C",x:0,y:0,z:0},{elem:"O",x:0,y:70,z:0},{elem:"O",x:60,y:-35,z:0},{elem:"O",x:-60,y:-35,z:0},{elem:"Ca",x:0,y:0,z:90,r:20,lpCount:0}], bonds: [[0,1,"double"],[0,2,"single"],[0,3,"single"]] },
    "MgCO3|碳酸鎂": { pg: "D3h", mp: "350 (分解)", bp: "-", desc: "<strong>碳酸鎂</strong><br>Mg²⁺ 位於碳酸根平面上方。", atoms: [{elem:"C",x:0,y:0,z:0},{elem:"O",x:0,y:70,z:0},{elem:"O",x:60,y:-35,z:0},{elem:"O",x:-60,y:-35,z:0},{elem:"Mg",x:0,y:0,z:90,r:18,lpCount:0}], bonds: [[0,1,"double"],[0,2,"single"],[0,3,"single"]] },
    "Na2CO3|碳酸鈉|蘇打": { pg: "D3h", mp: "851", bp: "-", desc: "<strong>碳酸鈉 (蘇打)</strong><br>兩個 Na⁺ 位於外側。", atoms: [{elem:"C",x:0,y:0,z:0},{elem:"O",x:0,y:70,z:0},{elem:"O",x:60,y:-35,z:0},{elem:"O",x:-60,y:-35,z:0},{elem:"Na",x:100,y:-20,z:0,r:15},{elem:"Na",x:-100,y:-20,z:0,r:15}], bonds: [[0,1,"double"],[0,2,"single"],[0,3,"single"]] },
    "K2CO3|碳酸鉀|草木灰": { pg: "D3h", mp: "891", bp: "-", desc: "<strong>碳酸鉀</strong><br>兩個 K⁺ 位於外側。", atoms: [{elem:"C",x:0,y:0,z:0},{elem:"O",x:0,y:70,z:0},{elem:"O",x:60,y:-35,z:0},{elem:"O",x:-60,y:-35,z:0},{elem:"K",x:100,y:-20,z:0,r:22},{elem:"K",x:-100,y:-20,z:0,r:22}], bonds: [[0,1,"double"],[0,2,"single"],[0,3,"single"]] },
    "NaHCO3|碳酸氫鈉|小蘇打": { pg: "Cs", mp: "50 (分解)", bp: "-", desc: "<strong>碳酸氫鈉</strong><br>Na⁺ 位於外側。", atoms: [{elem:"C",x:0,y:0,z:0},{elem:"O",x:0,y:70,z:0},{elem:"O",x:60,y:-35,z:0},{elem:"O",x:-60,y:-35,z:0},{elem:"H",x:-90,y:-10,z:0},{elem:"Na",x:100,y:-20,z:0,r:15}], bonds: [[0,1,"double"],[0,2],[0,3],[3,4]] }
});

addMol("HNO3|硝酸系列", "N", "sp²", ["平面三角形","Trigonal Planar"], "120°", "-42", "83", [], [], {
    "HNO3|硝酸": { pg: "Cs", mp: "-42", bp: "83", desc: "<strong>硝酸</strong><br>強酸及強氧化劑。光照易分解產生紅棕色 NO₂。", atoms: [{elem:"N",x:0,y:0,z:0,lpCount:0}, {elem:"O",x:0,y:68,z:0}, {elem:"O",x:-59,y:-34,z:0}, {elem:"O",x:59,y:-34,z:0,lpCount:2}, {elem:"H",x:90,y:-15,z:0}], bonds: [[0,1,"double"], [0,2,"coordinate"], [0,3], [3,4]] },
    "NO3-|硝酸根": { pg: "D3h", mp: "-", bp: "-", desc: "<strong>硝酸根</strong><br>具有高度對稱的平面結構 (共振)。", atoms: [{elem:"N",x:0,y:0,z:0,lpCount:0}, {elem:"O",x:0,y:68,z:0}, {elem:"O",x:-59,y:-34,z:0}, {elem:"O",x:59,y:-34,z:0}], bonds: [[0,1,"double"], [0,2,"coordinate"], [0,3]] },
    "KNO3|硝酸鉀|硝石": { pg: "D3h", mp: "334", bp: "400 (分解)", desc: "<strong>硝酸鉀</strong><br>俗稱硝石。K⁺ 位於結構上方。", atoms: [{elem:"N",x:0,y:0,z:0,lpCount:0}, {elem:"O",x:0,y:68,z:0}, {elem:"O",x:-59,y:-34,z:0}, {elem:"O",x:59,y:-34,z:0}, {elem:"K",x:0,y:0,z:90,r:22,lpCount:0}], bonds: [[0,1,"double"], [0,2,"coordinate"], [0,3]] },
    "NaNO3|硝酸鈉|智利硝石": { pg: "D3h", mp: "308", bp: "380 (分解)", desc: "<strong>硝酸鈉</strong><br>俗稱智利硝石。Na⁺ 位於結構上方。", atoms: [{elem:"N",x:0,y:0,z:0,lpCount:0}, {elem:"O",x:0,y:68,z:0}, {elem:"O",x:-59,y:-34,z:0}, {elem:"O",x:59,y:-34,z:0}, {elem:"Na",x:0,y:0,z:85,r:15,lpCount:0}], bonds: [[0,1,"double"], [0,2,"coordinate"], [0,3]] },
    "AgNO3|硝酸銀": { pg: "D3h", mp: "212", bp: "444 (分解)", desc: "<strong>硝酸銀</strong><br>Ag⁺ 位於結構上方。", atoms: [{elem:"N",x:0,y:0,z:0,lpCount:0}, {elem:"O",x:0,y:68,z:0}, {elem:"O",x:-59,y:-34,z:0}, {elem:"O",x:59,y:-34,z:0}, {elem:"Ag",x:0,y:0,z:90,r:18,lpCount:0}], bonds: [[0,1,"double"], [0,2,"coordinate"], [0,3]] },
    "Cu(NO3)2|硝酸銅": { pg: "D3h", mp: "114", bp: "170 (分解)", desc: "<strong>硝酸銅</strong><br>藍色晶體。Cu²⁺ 。", atoms: [{elem:"Cu",x:0,y:0,z:0,r:18,lpCount:0}, {elem:"N",x:-90,y:0,z:0,lpCount:0}, {elem:"O",x:-145,y:0,z:0}, {elem:"O",x:-60,y:45,z:35}, {elem:"O",x:-60,y:-45,z:-35}, {elem:"N",x:90,y:0,z:0,lpCount:0}, {elem:"O",x:145,y:0,z:0}, {elem:"O",x:60,y:45,z:35}, {elem:"O",x:60,y:-45,z:-35}], bonds: [[1,2,"double"],[1,3,"coordinate"],[1,4,"single"], [5,6,"double"],[5,7,"coordinate"],[5,8,"single"]] }
});

addMol("HNO2|亞硝酸系列", "N", "sp²", ["角形","Bent"], "111°", "-", "不穩定", [], [], {
    "HNO2|亞硝酸": { pg: "Cs", mp: "-", bp: "不穩定", desc: "<strong>亞硝酸</strong><br>弱酸，N原子上有一對孤對電子。", atoms: [{elem:"N",x:0,y:0,z:0,lpCount:1},{elem:"O",x:0,y:65,z:0},{elem:"O",x:60,y:-35,z:0},{elem:"H",x:90,y:-10,z:0}], bonds: [[0,1,"double"],[0,2],[2,3]] },
    "NO2-|亞硝酸根": { pg: "C2v", mp: "-", bp: "-", desc: "<strong>亞硝酸根</strong><br>常見的防腐劑成分(亞硝酸鹽)，結構呈V型。", atoms: [{elem:"N",x:0,y:0,z:0,lpCount:1},{elem:"O",x:0,y:65,z:0},{elem:"O",x:60,y:-35,z:0}], bonds: [[0,1,"double"],[0,2]] },
    "NaNO2|亞硝酸鈉": { pg: "C2v", mp: "271", bp: "320 (分解)", desc: "<strong>亞硝酸鈉</strong><br>Na⁺ 位於外側。", atoms: [{elem:"N",x:0,y:0,z:0,lpCount:1},{elem:"O",x:0,y:65,z:0},{elem:"O",x:60,y:-35,z:0},{elem:"Na",x:-80,y:0,z:0,r:15}], bonds: [[0,1,"double"],[0,2]] },
    "KNO2|亞硝酸鉀": { pg: "C2v", mp: "440 (分解)", bp: "-", desc: "<strong>亞硝酸鉀</strong><br>K⁺ 位於外側。", atoms: [{elem:"N",x:0,y:0,z:0,lpCount:1},{elem:"O",x:0,y:65,z:0},{elem:"O",x:60,y:-35,z:0},{elem:"K",x:-85,y:0,z:0,r:22}], bonds: [[0,1,"double"],[0,2]] }
});

addMol("H3PO4|磷酸系列", "P", "sp³", ["四面體","Tetrahedral"], "109.5°", "42.4", "213 (分解)", [], [], {
    "H3PO4|磷酸": { pg: "Cs", mp: "42.4", bp: "213 (分解)", desc: "<strong>磷酸</strong><br>三質子酸，含一個 P=O 與三個 P-OH。", atoms: [{elem:"P",x:0,y:0,z:0,lpCount:0},{elem:"O",x:0,y:68,z:0},{elem:"O",x:55,y:-30,z:35,lpCount:2},{elem:"O",x:-55,y:-30,z:35,lpCount:2},{elem:"O",x:0,y:-30,z:-60,lpCount:2},{elem:"H",x:80,y:-10,z:55},{elem:"H",x:-80,y:-10,z:55},{elem:"H",x:0,y:-10,z:-90}], bonds: [[0,1,"double"],[0,2],[0,3],[0,4],[2,5],[3,6],[4,7]] },
    "H2PO4-|磷酸二氫根": { pg: "C2v", mp: "-", bp: "-", desc: "<strong>磷酸二氫根</strong><br>帶 -1 價電荷。", atoms: [{elem:"P",x:0,y:0,z:0,lpCount:0},{elem:"O",x:0,y:68,z:0},{elem:"O",x:55,y:-30,z:35,lpCount:2},{elem:"O",x:-55,y:-30,z:35,lpCount:2},{elem:"O",x:0,y:-30,z:-60},{elem:"H",x:80,y:-10,z:55},{elem:"H",x:-80,y:-10,z:55}], bonds: [[0,1,"double"],[0,2],[0,3],[0,4],[2,5],[3,6]] },
    "HPO42-|磷酸氫根": { pg: "C3v", mp: "-", bp: "-", desc: "<strong>磷酸氫根</strong><br>帶 -2 價電荷。", atoms: [{elem:"P",x:0,y:0,z:0,lpCount:0},{elem:"O",x:0,y:68,z:0},{elem:"O",x:55,y:-30,z:35,lpCount:2},{elem:"O",x:-55,y:-30,z:35},{elem:"O",x:0,y:-30,z:-60},{elem:"H",x:80,y:-10,z:55}], bonds: [[0,1,"double"],[0,2],[0,3],[0,4],[2,5]] },
    "PO43-|磷酸根": { pg: "Td", mp: "-", bp: "-", desc: "<strong>磷酸根</strong><br>正四面體結構，四個 P-O 鍵長均等。", atoms: [{elem:"P",x:0,y:0,z:0,lpCount:0},{elem:"O",x:0,y:68,z:0},{elem:"O",x:55,y:-30,z:35},{elem:"O",x:-55,y:-30,z:35},{elem:"O",x:0,y:-30,z:-60}], bonds: [[0,1,"double"],[0,2],[0,3],[0,4]] },
    "Ca3(PO4)2|磷酸鈣": { pg: "Td", mp: "1670", bp: "-", desc: "<strong>磷酸鈣</strong><br>難溶於水，變量原料。", atoms: [{elem:"P",x:0,y:0,z:0,lpCount:0},{elem:"O",x:0,y:68,z:0},{elem:"O",x:55,y:-30,z:35},{elem:"O",x:-55,y:-30,z:35},{elem:"O",x:0,y:-30,z:-60},{elem:"Ca",x:100,y:40,z:0,r:20},{elem:"Ca",x:-100,y:40,z:0,r:20},{elem:"Ca",x:0,y:-100,z:0,r:20}], bonds: [[0,1,"double"],[0,2],[0,3],[0,4]] },
    "Na3PO4|磷酸鈉": { pg: "Td", mp: "1583", bp: "-", desc: "<strong>磷酸鈉</strong><br>強鹼性鹽類。", atoms: [{elem:"P",x:0,y:0,z:0,lpCount:0},{elem:"O",x:0,y:68,z:0},{elem:"O",x:55,y:-30,z:35},{elem:"O",x:-55,y:-30,z:35},{elem:"O",x:0,y:-30,z:-60},{elem:"Na",x:90,y:30,z:0,r:15},{elem:"Na",x:-90,y:30,z:0,r:15},{elem:"Na",x:0,y:-90,z:0,r:15}], bonds: [[0,1,"double"],[0,2],[0,3],[0,4]] },
    "Ca(H2PO4)2|磷酸二氫鈣": { pg: "C2v", mp: "109 (分解)", bp: "-", desc: "<strong>磷酸二氫鈣</strong><br>肥料成分。", atoms: [{elem:"Ca",x:0,y:0,z:0,r:20}, {elem:"P",x:-100,y:0,z:0}, {elem:"O",x:-100,y:65,z:0}, {elem:"O",x:-100,y:-30,z:55}, {elem:"O",x:-145,y:-30,z:-30}, {elem:"O",x:-55,y:-30,z:-30}, {elem:"H",x:-145,y:-60,z:55}, {elem:"H",x:-175,y:-10,z:-30}, {elem:"P",x:100,y:0,z:0}, {elem:"O",x:100,y:65,z:0}, {elem:"O",x:100,y:-30,z:55}, {elem:"O",x:145,y:-30,z:-30}, {elem:"O",x:55,y:-30,z:-30}, {elem:"H",x:145,y:-60,z:55}, {elem:"H",x:175,y:-10,z:-30}], bonds: [[1,2,"double"],[1,3],[1,4],[1,5],[3,6],[4,7], [8,9,"double"],[8,10],[8,11],[8,12],[10,13],[11,14]] }
});

addMol("H3PO3|亞磷酸系列", "P", "sp³", ["四面體","Tetrahedral"], "109.5°", "73.6", "200 (分解)", [], [], {
    "H3PO3|亞磷酸": { pg: "Cs", mp: "73.6", bp: "200 (分解)", desc: "<strong>亞磷酸</strong><br>二質子酸，含一個 P-H 鍵 (不解離) 與兩個 P-OH。", atoms: [{elem:"P",x:0,y:0,z:0,lpCount:0},{elem:"O",x:0,y:68,z:0},{elem:"O",x:55,y:-30,z:35,lpCount:2},{elem:"O",x:-55,y:-30,z:35,lpCount:2},{elem:"H",x:0,y:-40,z:-60},{elem:"H",x:90,y:-10,z:60},{elem:"H",x:-90,y:-10,z:60}], bonds: [[0,1,"double"],[0,2],[0,3],[0,4],[2,5],[3,6]] },
    "H2PO3-|亞磷酸氫根": { pg: "Cs", mp: "-", bp: "-", desc: "<strong>亞磷酸二氫根</strong><br>帶 -1 價電荷，P-H 鍵保留。", atoms: [{elem:"P",x:0,y:0,z:0,lpCount:0},{elem:"O",x:0,y:68,z:0},{elem:"O",x:55,y:-30,z:35,lpCount:2},{elem:"O",x:-55,y:-30,z:35,lpCount:2},{elem:"H",x:0,y:-40,z:-60},{elem:"H",x:90,y:-10,z:60}], bonds: [[0,1,"double"],[0,2],[0,3],[0,4],[2,5]] },
    "HPO32-|亞磷酸根": { pg: "C3v", mp: "-", bp: "-", desc: "<strong>亞磷酸氫根 (亞磷酸根)</strong><br>帶 -2 價電荷，P-H 鍵通常不解離。", atoms: [{elem:"P",x:0,y:0,z:0,lpCount:0},{elem:"O",x:0,y:68,z:0},{elem:"O",x:55,y:-30,z:35},{elem:"O",x:-55,y:-30,z:35},{elem:"H",x:0,y:-40,z:-60}], bonds: [[0,1,"double"],[0,2],[0,3],[0,4]] },
    "Na2HPO3|亞磷酸鈉": { pg: "C3v", mp: "-", bp: "-", desc: "<strong>亞磷酸鈉</strong><br>正鹽，P 直接連有一個 H。", atoms: [{elem:"P",x:0,y:0,z:0,lpCount:0},{elem:"O",x:0,y:68,z:0},{elem:"O",x:55,y:-30,z:35},{elem:"O",x:-55,y:-30,z:35},{elem:"H",x:0,y:-40,z:-60},{elem:"Na",x:90,y:20,z:0,r:15},{elem:"Na",x:-90,y:20,z:0,r:15}], bonds: [[0,1,"double"],[0,2],[0,3],[0,4]] }
});

addMol("H3PO2|次磷酸系列", "P", "sp³", ["四面體","Tetrahedral"], "109.5°", "26.5", "130 (分解)", [], [], {
    "H3PO2|次磷酸": { pg: "Cs", mp: "26.5", bp: "130 (分解)", desc: "<strong>次磷酸</strong><br>單質子酸，含兩個 P-H 鍵與一個 P-OH。", atoms: [{elem:"P",x:0,y:0,z:0,lpCount:0},{elem:"O",x:0,y:68,z:0},{elem:"O",x:0,y:-30,z:-60,lpCount:2},{elem:"H",x:55,y:-35,z:35},{elem:"H",x:-55,y:-35,z:35},{elem:"H",x:0,y:-10,z:-100}], bonds: [[0,1,"double"],[0,2],[0,3],[0,4],[2,5]] },
    "H2PO2-|次磷酸根": { pg: "C2v", mp: "-", bp: "-", desc: "<strong>次磷酸根</strong><br>帶 -1 價電荷。", atoms: [{elem:"P",x:0,y:0,z:0,lpCount:0},{elem:"O",x:0,y:68,z:0},{elem:"O",x:0,y:-30,z:-60},{elem:"H",x:55,y:-35,z:35},{elem:"H",x:-55,y:-35,z:35}], bonds: [[0,1,"double"],[0,2],[0,3],[0,4]] },
    "NaH2PO2|次磷酸鈉": { pg: "C2v", mp: "90 (一水合)", bp: "-", desc: "<strong>次磷酸鈉</strong><br>強還原劑。", atoms: [{elem:"P",x:0,y:0,z:0,lpCount:0},{elem:"O",x:0,y:68,z:0},{elem:"O",x:0,y:-30,z:-60},{elem:"H",x:55,y:-35,z:35},{elem:"H",x:-55,y:-35,z:35},{elem:"Na",x:-85,y:0,z:0,r:15}], bonds: [[0,1,"double"],[0,2],[0,3],[0,4]] }
});

addMol("HClO4|過氯酸系列", "Cl", "sp³", ["四面體","Tetrahedral"], "109.5°", "-112", "19 (分解)", [], [], {
    "HClO4|過氯酸": { pg: "Cs", mp: "-112", bp: "19 (分解)", desc: "<strong>過氯酸</strong><br>最強無機酸之一，正四面體結構。氯原子與三個氧形成雙鍵，與一個羥基形成單鍵。", atoms: [{elem:"Cl",x:0,y:0,z:0,lpCount:0},{elem:"O",x:0,y:68,z:0},{elem:"O",x:58,y:-25,z:35,lpCount:2},{elem:"O",x:-58,y:-25,z:35,lpCount:2},{elem:"O",x:0,y:-25,z:-65,lpCount:2},{elem:"H",x:0,y:-5,z:-105}], bonds: [[0,1,"double"],[0,2,"double"],[0,3,"double"],[0,4,"single"],[4,5]] },
    "ClO4-|過氯酸根": { pg: "Td", mp: "-", bp: "-", desc: "<strong>過氯酸根</strong><br>化學性質穩定，四個 Cl-O 鍵長因共振而均等 (-1價)。", atoms: [{elem:"Cl",x:0,y:0,z:0,lpCount:0},{elem:"O",x:0,y:68,z:0},{elem:"O",x:58,y:-25,z:35},{elem:"O",x:-58,y:-25,z:35},{elem:"O",x:0,y:-25,z:-65}], bonds: [[0,1,"double"],[0,2,"double"],[0,3,"double"],[0,4,"single"]] },
    "Mg(ClO4)2|過氯酸鎂": { pg: "Td", mp: "251", bp: "-", desc: "<strong>過氯酸鎂</strong><br>極強的脫水劑（乾燥劑）。", atoms: [{elem:"Mg",x:0,y:0,z:0,r:20,lpCount:0}, {elem:"Cl",x:-130,y:0,z:0,lpCount:0},{elem:"O",x:-130,y:68,z:0},{elem:"O",x:-72,y:-25,z:35},{elem:"O",x:-188,y:-25,z:35},{elem:"O",x:-130,y:-25,z:-65}, {elem:"Cl",x:130,y:0,z:0,lpCount:0},{elem:"O",x:130,y:68,z:0},{elem:"O",x:72,y:-25,z:35},{elem:"O",x:188,y:-25,z:35},{elem:"O",x:130,y:-25,z:-65}], bonds: [[1,2,"double"],[1,3,"double"],[1,4,"double"],[1,5,"single"], [6,7,"double"],[6,8,"double"],[6,9,"double"],[6,10,"single"]] },
    "KClO4|過氯酸鉀": { pg: "Td", mp: "610 (分解)", bp: "-", desc: "<strong>過氯酸鉀</strong><br>強氧化劑，用於煙火（紫色火焰）。", atoms: [{elem:"Cl",x:0,y:0,z:0,lpCount:0},{elem:"O",x:0,y:68,z:0},{elem:"O",x:58,y:-25,z:35},{elem:"O",x:-58,y:-25,z:35},{elem:"O",x:0,y:-25,z:-65},{elem:"K",x:0,y:0,z:95,r:22}], bonds: [[0,1,"double"],[0,2,"double"],[0,3,"double"],[0,4,"single"]] },
    "NH4ClO4|過氯酸銨": { pg: "Td", mp: "240 (分解)", bp: "-", desc: "<strong>過氯酸銨 (AP)</strong><br>固體火箭燃料氧化劑。", atoms: [{elem:"Cl",x:0,y:0,z:0,lpCount:0},{elem:"O",x:0,y:68,z:0},{elem:"O",x:58,y:-25,z:35},{elem:"O",x:-58,y:-25,z:35},{elem:"O",x:0,y:-25,z:-65},{elem:"N",x:110,y:0,z:0,r:18},{elem:"H",x:110,y:40,z:0},{elem:"H",x:110,y:-20,z:35},{elem:"H",x:110,y:-20,z:-35},{elem:"H",x:145,y:0,z:0}], bonds: [[0,1,"double"],[0,2,"double"],[0,3,"double"],[0,4,"single"],[5,6],[5,7],[5,8],[5,9]] }
});

addMol("HClO3|氯酸系列", "Cl", "sp³", ["角錐形","Pyramidal"], "107°", "-20", "分解", [], [], {
    "HClO3|氯酸": { pg: "Cs", mp: "-20", bp: "分解", desc: "<strong>氯酸</strong><br>強酸，具有強氧化性，中心有一對孤對電子。", atoms: [{elem:"Cl",x:0,y:15,z:0,lpCount:1},{elem:"O",x:0,y:-40,z:50},{elem:"O",x:48,y:-40,z:-28},{elem:"O",x:-48,y:-40,z:-28,lpCount:2},{elem:"H",x:-90,y:-20,z:-55}], bonds: [[0,1,"double"],[0,2,"double"],[0,3,"single"],[3,4]] },
    "ClO3-|氯酸根": { pg: "C3v", mp: "-", bp: "-", desc: "<strong>氯酸根</strong><br>三角錐形結構，常用於火藥與炸藥。", atoms: [{elem:"Cl",x:0,y:15,z:0,lpCount:1},{elem:"O",x:0,y:-40,z:50},{elem:"O",x:48,y:-40,z:-28},{elem:"O",x:-48,y:-40,z:-28}], bonds: [[0,1,"double"],[0,2,"double"],[0,3,"single"]] },
    "KClO3|氯酸鉀": { pg: "C3v", mp: "356", bp: "400 (分解)", desc: "<strong>氯酸鉀</strong><br>強氧化劑，受熱分解產生氧氣。", atoms: [{elem:"Cl",x:0,y:15,z:0,lpCount:1},{elem:"O",x:0,y:-40,z:50},{elem:"O",x:48,y:-40,z:-28},{elem:"O",x:-48,y:-40,z:-28},{elem:"K",x:0,y:0,z:85,r:22}], bonds: [[0,1,"double"],[0,2,"double"],[0,3,"single"]] },
    "NaClO3|氯酸鈉": { pg: "C3v", mp: "248", bp: "300 (分解)", desc: "<strong>氯酸鈉</strong><br>工業漂白與除草劑原料。", atoms: [{elem:"Cl",x:0,y:15,z:0,lpCount:1},{elem:"O",x:0,y:-40,z:50},{elem:"O",x:48,y:-40,z:-28},{elem:"O",x:-48,y:-40,z:-28},{elem:"Na",x:0,y:0,z:80,r:15}], bonds: [[0,1,"double"],[0,2,"double"],[0,3,"single"]] }
});

addMol("HClO2|亞氯酸系列", "Cl", "sp³", ["角形","Bent"], "111°", "-", "不穩定", [], [], {
    "HClO2|亞氯酸": { pg: "Cs", mp: "-", bp: "不穩定", desc: "<strong>亞氯酸</strong><br>弱酸，結構呈V型，中心有兩對孤對電子。", atoms: [{elem:"Cl",x:0,y:5,z:0,lpCount:2},{elem:"O",x:55,y:-35,z:0},{elem:"O",x:-55,y:-35,z:0,lpCount:2},{elem:"H",x:-90,y:-20,z:0}], bonds: [[0,1,"double"],[0,2,"single"],[2,3]] },
    "ClO2-|亞氯酸根": { pg: "C2v", mp: "-", bp: "-", desc: "<strong>亞氯酸根</strong><br>V型結構，常用於漂白劑。", atoms: [{elem:"Cl",x:0,y:5,z:0,lpCount:2},{elem:"O",x:55,y:-35,z:0},{elem:"O",x:-55,y:-35,z:0}], bonds: [[0,1,"double"],[0,2,"single"]] },
    "NaClO2|亞氯酸鈉": { pg: "C2v", mp: "170 (分解)", bp: "-", desc: "<strong>亞氯酸鈉</strong><br>高效漂白劑，反應可生成二氧化氯 (ClO₂)。", atoms: [{elem:"Cl",x:0,y:5,z:0,lpCount:2},{elem:"O",x:55,y:-35,z:0},{elem:"O",x:-55,y:-35,z:0},{elem:"Na",x:-90,y:0,z:0,r:15}], bonds: [[0,1,"double"],[0,2,"single"]] }
});

addMol("HClO|次氯酸系列", "O", "sp³", ["角形","Bent"], "104.5°", "-", "不穩定", [], [], {
    "HClO|次氯酸": { pg: "Cs", mp: "-", bp: "不穩定", desc: "<strong>次氯酸</strong><br>弱酸，殺菌力強，結構 H-O-Cl。", atoms: [{elem:"O",x:0,y:10,z:0,lpCount:2},{elem:"Cl",x:65,y:-25,z:0},{elem:"H",x:-35,y:-20,z:0}], bonds: [[0,1],[0,2]] },
    "ClO-|次氯酸根": { pg: "Cinfv", mp: "-", bp: "-", desc: "<strong>次氯酸根</strong><br>漂白水有效成分。", atoms: [{elem:"Cl",x:-35,y:0,z:0,lpCount:3},{elem:"O",x:35,y:0,z:0,lpCount:3}], bonds: [[0,1]] },
    "NaClO|次氯酸鈉|漂白水": { pg: "Cinfv", mp: "18 (五水合)", bp: "分解", desc: "<strong>次氯酸鈉 (漂白水)</strong><br>家用漂白劑。Na⁺ 與 ClO⁻ 之間為離子鍵。", atoms: [{elem:"Cl",x:-35,y:0,z:0,lpCount:3},{elem:"O",x:35,y:0,z:0,lpCount:3},{elem:"Na",x:85,y:0,z:0,r:15}], bonds: [[0,1]] },
    "Ca(ClO)2|次氯酸鈣|漂白粉": { pg: "Cinfv", mp: "100 (分解)", bp: "-", desc: "<strong>次氯酸鈣</strong><br>漂白粉主要成分。", atoms: [{elem:"Cl",x:-55,y:0,z:0,lpCount:3},{elem:"O",x:15,y:0,z:0,lpCount:3},{elem:"Ca",x:60,y:0,z:0,r:20},{elem:"O",x:105,y:0,z:0,lpCount:3},{elem:"Cl",x:175,y:0,z:0,lpCount:3}], bonds: [[0,1], [3,4]] }
});

addMol("HBrO3|溴酸系列", "Br", "sp³", ["角錐形","Pyramidal"], "107°", "-", "不穩定", [], [], {
    "HBrO3|溴酸": { pg: "Cs", mp: "-", bp: "不穩定", desc: "<strong>溴酸</strong><br>強酸，中心有一對孤對電子。", atoms: [{elem:"Br",x:0,y:15,z:0,lpCount:1},{elem:"O",x:0,y:-40,z:50},{elem:"O",x:48,y:-40,z:-28},{elem:"O",x:-48,y:-40,z:-28,lpCount:2},{elem:"H",x:-90,y:-20,z:-55}], bonds: [[0,1,"double"],[0,2,"double"],[0,3,"single"],[3,4]] },
    "BrO3-|溴酸根": { pg: "C3v", mp: "-", bp: "-", desc: "<strong>溴酸根</strong><br>三角錐形結構。", atoms: [{elem:"Br",x:0,y:15,z:0,lpCount:1},{elem:"O",x:0,y:-40,z:50},{elem:"O",x:48,y:-40,z:-28},{elem:"O",x:-48,y:-40,z:-28}], bonds: [[0,1,"double"],[0,2,"double"],[0,3,"single"]] },
    "KBrO3|溴酸鉀": { pg: "C3v", mp: "350 (分解)", bp: "-", desc: "<strong>溴酸鉀</strong><br>強氧化劑，K⁺ 位於外側。", atoms: [{elem:"Br",x:0,y:15,z:0,lpCount:1},{elem:"O",x:0,y:-40,z:50},{elem:"O",x:48,y:-40,z:-28},{elem:"O",x:-48,y:-40,z:-28},{elem:"K",x:0,y:60,z:0,r:22}], bonds: [[0,1,"double"],[0,2,"double"],[0,3,"single"]] },
    "AgBrO3|溴酸銀": { pg: "C3v", mp: "-", bp: "-", desc: "<strong>溴酸銀</strong><br>難溶於水的白色固體。", atoms: [{elem:"Br",x:0,y:15,z:0,lpCount:1},{elem:"O",x:0,y:-40,z:50},{elem:"O",x:48,y:-40,z:-28},{elem:"O",x:-48,y:-40,z:-28},{elem:"Ag",x:0,y:60,z:0,r:18}], bonds: [[0,1,"double"],[0,2,"double"],[0,3,"single"]] }
});

addMol("HIO3|碘酸系列", "I", "sp³", ["角錐形","Pyramidal"], "107°", "110", "分解", [], [], {
    "HIO3|碘酸": { pg: "Cs", mp: "110", bp: "分解", desc: "<strong>碘酸</strong><br>穩定的白色固體，強酸。", atoms: [{elem:"I",x:0,y:15,z:0,lpCount:1},{elem:"O",x:0,y:-40,z:50},{elem:"O",x:48,y:-40,z:-28},{elem:"O",x:-48,y:-40,z:-28,lpCount:2},{elem:"H",x:-90,y:-20,z:-55}], bonds: [[0,1,"double"],[0,2,"double"],[0,3,"single"],[3,4]] },
    "IO3-|碘酸根": { pg: "C3v", mp: "-", bp: "-", desc: "<strong>碘酸根</strong><br>三角錐形結構。", atoms: [{elem:"I",x:0,y:15,z:0,lpCount:1},{elem:"O",x:0,y:-40,z:50},{elem:"O",x:48,y:-40,z:-28},{elem:"O",x:-48,y:-40,z:-28}], bonds: [[0,1,"double"],[0,2,"double"],[0,3,"single"]] },
    "KIO3|碘酸鉀": { pg: "C3v", mp: "560 (分解)", bp: "-", desc: "<strong>碘酸鉀</strong><br>食鹽加碘成分，K⁺ 位於外側。", atoms: [{elem:"I",x:0,y:15,z:0,lpCount:1},{elem:"O",x:0,y:-40,z:50},{elem:"O",x:48,y:-40,z:-28},{elem:"O",x:-48,y:-40,z:-28},{elem:"K",x:0,y:60,z:0,r:22}], bonds: [[0,1,"double"],[0,2,"double"],[0,3,"single"]] }
});

// --- 11. 簡單有機分子與衍生物 (鍵長修正: C-H=50, C-C=70, C-N=70, C=O=68, C-Cl=75) ---
addMol("CH3NO2|硝基甲烷", "C", "sp3", ["正四面體","Tetrahedral"], "109.5°", "-29", "101.2", [{elem:"O",x:-65,y:50,z:0,lpCount:3},{elem:"O",x:-65,y:-50,z:0},{elem:"N",x:-35,y:0,z:0,lpCount:0},{elem:"C",x:35,y:0,z:0},{elem:"H",x:55,y:-35,z:-25},{elem:"H",x:55,y:35,z:-25},{elem:"H",x:55,y:0,z:45}], [[2,0,"coordinate"],[1,2,"double"],[2,3],[3,4],[3,5],[3,6]], null, null, "Cs");
addMol("C2H6|乙烷", "C", "sp³", ["四面體連結","Tetrahedral"], "109.5°", "-182.8", "-88.6", [{elem:"C",x:-35,y:0,z:0}, {elem:"C",x:35,y:0,z:0}, {elem:"H",x:-65,y:35,z:0}, {elem:"H",x:-65,y:-25,z:30}, {elem:"H",x:-65,y:-25,z:-30}, {elem:"H",x:65,y:-35,z:0}, {elem:"H",x:65,y:25,z:30}, {elem:"H",x:65,y:25,z:-30}], [[0,1], [0,2], [0,3], [0,4], [1,5], [1,6], [1,7]], null, null, "D3d");
addMol("CH3Cl|一氯甲烷|氯甲烷", "C", "sp3", ["正四面體","Tetrahedral"], "109.5°", "-97.4", "-24.2", [{elem:"Cl",x:75,y:0,z:0},{elem:"C",x:0,y:0,z:0},{elem:"H",x:-30,y:30,z:30},{elem:"H",x:-30,y:-40,z:10},{elem:"H",x:-30,y:10,z:-40}], [[0,1],[1,2],[1,3],[1,4]], null, null, "C3v");
addMol("C2H4|乙烯", "C", "sp²", ["平面","Planar"], "120°", "-169.2", "-103.7", [{elem:"C",x:-32,y:0,z:0}, {elem:"C",x:32,y:0,z:0}, {elem:"H",x:-67,y:45,z:0}, {elem:"H",x:-67,y:-45,z:0}, {elem:"H",x:67,y:45,z:0}, {elem:"H",x:67,y:-45,z:0}], [[0,1,"double"], [0,2], [0,3], [1,4], [1,5]], null, null, "D2h");
addMol("C2H2|乙炔", "C", "sp", ["直線型","Linear"], "180°", "-80.8", "-84 (昇華)", [{elem:"C",x:-30,y:0,z:0}, {elem:"C",x:30,y:0,z:0}, {elem:"H",x:-80,y:0,z:0}, {elem:"H",x:80,y:0,z:0}], [[0,1,"triple"], [0,2], [1,3]], null, null, "Dinfh");

addMol("C2H2", "C", "sp", ["直線型","Linear"], "180°", "-80.8", "-84", [], [], {
    "C2H2|乙炔": { pg: "Dinfh", mp: "-80.8", bp: "-84", atoms: [{elem:"C",x:-30,y:0,z:0},{elem:"C",x:30,y:0,z:0},{elem:"H",x:-80,y:0,z:0},{elem:"H",x:80,y:0,z:0}], bonds: [[0,1,"triple"], [0,2], [1,3]] }
});

addMol("N2H4|聯氨|肼", "N", "sp³", ["扭轉型","Gauche"], "107°", "2", "114", [{elem:"N",x:-35,y:0,z:0}, {elem:"N",x:35,y:0,z:0}, {elem:"H",x:-60,y:35,z:25}, {elem:"H",x:-60,y:-35,z:25}, {elem:"H",x:60,y:35,z:-25}, {elem:"H",x:60,y:-35,z:-25}], [[0,1], [0,2], [0,3], [1,4], [1,5]], null, null, "C2");

addMol("C3H6O|丙醛/丙酮 (同分異構)", "C", "sp²", ["平面/四面體","Mixed"], "120°", "-81", "48", [], [], 
    {
        "C2H5CHO|丙醛": { pg: "Cs", mp: "-81", bp: "48", atoms: [{elem:"O",x:-86,y:35,z:0},{elem:"C",x:5,y:-25,z:0},{elem:"C",x:45,y:35,z:0},{elem:"C",x:-60,y:-15,z:0},{elem:"H",x:15,y:-55,z:40},{elem:"H",x:15,y:-55,z:-40},{elem:"H",x:35,y:65,z:-40},{elem:"H",x:35,y:65,z:40},{elem:"H",x:95,y:25,z:0},{elem:"H",x:-85,y:-60,z:0}],  bonds: [[0,3,"double"],[1,2],[1,3],[1,4],[1,5],[2,6],[2,7],[2,8],[3,9]] },
        "CH3COCH3|丙酮": { pg: "C2v", mp: "-94.7", bp: "56.1", atoms: [{elem:"O",x:0,y:85,z:0},{elem:"C",x:0,y:20,z:0},{elem:"C",x:60,y:-15,z:0},{elem:"C",x:-60,y:-15,z:0},{elem:"H",x:60,y:-45,z:-40},{elem:"H",x:60,y:-45,z:40},{elem:"H",x:100,y:15,z:0},{elem:"H",x:-100,y:15,z:0},{elem:"H",x:-60,y:-45,z:-40},{elem:"H",x:-60,y:-45,z:40}], bonds: [[0,1,"double"],[1,2],[1,3],[2,4],[2,5],[2,6],[3,7],[3,8],[3,9]] }
    }
);

addMol("CH3NH2|甲胺|甲基胺", "N", "sp³", ["角錐形 (N端)","Pyramidal"], "107°", "-93", "-6.3", [{elem:"N",x:0,y:15,z:0,lp3d:[{x:0,y:1,z:0}]}, {elem:"H",x:30,y:-25,z:20}, {elem:"H",x:-30,y:-25,z:20}, {elem:"C",x:0,y:-25,z:-35}, {elem:"H",x:0,y:5,z:-65}, {elem:"H",x:25,y:-40,z:-55}, {elem:"H",x:-25,y:-40,z:-55}], [[0,1],[0,2],[0,3],[3,4],[3,5],[3,6]], null, null, "Cs");
addMol("CH3OCH3|甲醚|二甲醚", "O", "sp³", ["角形 (中心)","Bent"], "111°", "-141", "-24.8", [{elem:"O",x:0,y:15,z:0,lp3d:[{x:0,y:1,z:1},{x:0,y:1,z:-1}]}, {elem:"C",x:50,y:-20,z:0}, {elem:"C",x:-50,y:-20,z:0}, {elem:"H",x:80,y:0,z:0}, {elem:"H",x:50,y:-50,z:25}, {elem:"H",x:50,y:-50,z:-25}, {elem:"H",x:-80,y:0,z:0}, {elem:"H",x:-50,y:-50,z:25}, {elem:"H",x:-50,y:-50,z:-25}], [[0,1],[0,2],[1,3],[1,4],[1,5],[2,6],[2,7],[2,8]], null, null, "C2v");
addMol("HCOOH|甲酸|蟻酸", "C", "sp²", ["平面","Planar"], "120°", "8.4", "100.8", [{elem:"C",x:0,y:0,z:0}, {elem:"O",x:0,y:60,z:0}, {elem:"O",x:50,y:-35,z:0}, {elem:"H",x:-50,y:-35,z:0}, {elem:"H",x:80,y:-15,z:0}], [[0,1,"double"],[0,2],[0,3],[2,4]], null, null, "Cs");
addMol("CH2Cl2|二氯甲烷", "C", "sp³", ["四面體","Tetrahedral"], "109.5°", "-97", "39.6", getTetra("C","H", 50).map((a,i)=>i===1||i===2?{...a,elem:"Cl", x:a.x*1.5, y:a.y*1.5, z:a.z*1.5}:a), [[0,1],[0,2],[0,3],[0,4]], null, null, "C2v");
addMol("CHCl3|三氯甲烷|氯仿", "C", "sp³", ["四面體","Tetrahedral"], "109.5°", "-63.5", "61.2", getTetra("C","Cl", 75).map((a,i)=>i===4?{...a,elem:"H", x:a.x*0.66, y:a.y*0.66, z:a.z*0.66}:a), [[0,1],[0,2],[0,3],[0,4]], null, null, "C3v");
addMol("HCN|氰化氫|氫氰酸", "C", "sp", ["直線型","Linear"], "180°", "-13.3", "26", [{elem:"C",x:0,y:0,z:0}, {elem:"N",x:60,y:0,z:0,lpCount:1}, {elem:"H",x:-50,y:0,z:0}], [[0,1,"triple"],[0,2]], null, null, "Cinfv");
addMol("CH3CN|乙腈|氰甲烷", "C", "sp", ["直線型 (CN端)","Linear"], "180°", "-45", "82", [{elem:"C",x:0,y:0,z:0}, {elem:"N",x:60,y:0,z:0,lpCount:1}, {elem:"C",x:-70,y:0,z:0}, {elem:"H",x:-100,y:25,z:0}, {elem:"H",x:-100,y:-15,z:20}, {elem:"H",x:-100,y:-15,z:-20}], [[0,1,"triple"],[0,2],[2,3],[2,4],[2,5]], null, null, "C3v");
addMol("CO(NH2)2|尿素", "C", "sp²", ["平面","Planar"], "120°", "132.7", "分解", [{elem:"C",x:0,y:0,z:0}, {elem:"O",x:0,y:60,z:0}, {elem:"N",x:-50,y:-35,z:0}, {elem:"N",x:50,y:-35,z:0}, {elem:"H",x:-80,y:-15,z:0}, {elem:"H",x:-50,y:-70,z:0}, {elem:"H",x:80,y:-15,z:0}, {elem:"H",x:50,y:-70,z:0}], [[0,1,"double"],[0,2],[0,3],[2,4],[2,5],[3,6],[3,7]], null, null, "C2v");
addMol("SOCl2|亞硫醯氯|二氯亞碸", "S", "sp³", ["角錐形","Pyramidal"], "106°", "-104.5", "76", [{elem:"S",x:0,y:15,z:0,lp3d:[{x:0,y:1,z:0}]}, {elem:"O",x:0,y:-35,z:50}, {elem:"Cl",x:60,y:-35,z:-30}, {elem:"Cl",x:-60,y:-35,z:-30}], [[0,1,"double"],[0,2],[0,3]], null, null, "Cs");
addMol("POCl3|三氯氧化磷|磷醯氯", "P", "sp³", ["四面體","Tetrahedral"], "109.5°", "1.25", "105.8", [{elem:"P",x:0,y:0,z:0}, {elem:"O",x:0,y:65,z:0}, {elem:"Cl",x:55,y:-35,z:35}, {elem:"Cl",x:-55,y:-35,z:35}, {elem:"Cl",x:0,y:0,z:-70}], [[0,1,"double"],[0,2],[0,3],[0,4]], null, null, "C3v");
// C60 (使用 IIFE 生成資料並添加名稱)
(function(){
    const phi=(1+Math.sqrt(5))/2,scale=28;let rawVerts=[];
    const groups=[[0,1,3*phi],[1,2+phi,2*phi],[phi,2,2*phi+1]];
    function addPermutations(u,v,w){
        const cycles=[[u,v,w],[v,w,u],[w,u,v]];
        cycles.forEach(p=>{
            for(let i=0;i<8;i++){
                let x=p[0]*((i&1)?1:-1),y=p[1]*((i&2)?1:-1),z=p[2]*((i&4)?1:-1);
                rawVerts.push({x,y,z});
            }
        });
    }
    groups.forEach(g=>addPermutations(g[0],g[1],g[2]));
    const atoms=[],threshold=0.1;
    rawVerts.forEach(v=>{
        if(!atoms.some(a=>Math.abs(a.x-v.x*scale)<threshold&&Math.abs(a.y-v.y*scale)<threshold&&Math.abs(a.z-v.z*scale)<threshold)){
            atoms.push({elem:"C",x:v.x*scale,y:v.y*scale,z:v.z*scale,lpCount:0,r:5});
        }
    });
    const bonds=[],bondSet=new Set();
    for(let i=0;i<atoms.length;i++){
        let distances=[];
        for(let j=0;j<atoms.length;j++){
            if(i===j)continue;
            let d=(atoms[i].x-atoms[j].x)**2+(atoms[i].y-atoms[j].y)**2+(atoms[i].z-atoms[j].z)**2;
            distances.push({id:j,dist:d});
        }
        distances.sort((a,b)=>a.dist-b.dist);
        for(let k=0;k<3;k++){
            let neighbor=distances[k].id,id1=Math.min(i,neighbor),id2=Math.max(i,neighbor),key=`${id1}-${id2}`;
            if(!bondSet.has(key)){
                let type=(k===0)?"double":"single";
                bonds.push([id1,id2,type]);
                bondSet.add(key);
            }
        }
    }
    // [修正] 加入中文名稱搜尋
    addMol("C60|碳60|碳六十|富勒烯|足球烯","C","sp²",["球狀 (12個五邊形, 20個六邊形)","Truncated Icosahedron"],"120°", "> 280 (昇華)", "-", atoms,bonds);
})();

// --- 12. 有同分異構物的有機分子與其他有機化合物 (全面修正：Key 統一為 分子式|中文名稱) ---
// N2F2 (二氟二氮)
addMol("N2F2", "N", "sp²", ["平面","Planar"], "120°", "-165", "-105", [], [], {
    "N2F2|順式-二氟二氮": { pg: "C2v", mp: "-165", bp: "-105", atoms: [{elem:"N",x:-32,y:0,z:0},{elem:"N",x:32,y:0,z:0},{elem:"F",x:-67,y:60,z:0},{elem:"F",x:67,y:60,z:0}], bonds: [[0,1,"double"],[0,2],[1,3]] },
    "N2F2|反式-二氟二氮": { pg: "C2h", mp: "-172", bp: "-111", atoms: [{elem:"N",x:-32,y:0,z:0},{elem:"N",x:32,y:0,z:0},{elem:"F",x:-67,y:60,z:0},{elem:"F",x:67,y:-60,z:0}], bonds: [[0,1,"double"],[0,2],[1,3]] }
});

// C2H2Cl2 (二氯乙烯)
addMol("C2H2Cl2|二氯乙烯", "C", "sp²", ["平面","Planar"], "120°", "-81.5", "60.3", [], [], {
    "C2H2Cl2|順-1,2-二氯乙烯": { pg: "C2v", mp: "-81.5", bp: "60.3", atoms: [{elem:"C",x:-32,y:0,z:0},{elem:"C",x:32,y:0,z:0},{elem:"Cl",x:-69,y:65,z:0},{elem:"Cl",x:69,y:65,z:0},{elem:"H",x:-57,y:-43,z:0},{elem:"H",x:57,y:-43,z:0}], bonds: [[0,1,"double"],[0,2],[0,4],[1,3],[1,5]] },
    "C2H2Cl2|反-1,2-二氯乙烯": { pg: "C2h", mp: "-49.4", bp: "47.5", atoms: [{elem:"C",x:-32,y:0,z:0},{elem:"C",x:32,y:0,z:0},{elem:"Cl",x:-69,y:65,z:0},{elem:"Cl",x:69,y:-65,z:0},{elem:"H",x:-57,y:-43,z:0},{elem:"H",x:57,y:43,z:0}], bonds: [[0,1,"double"],[0,2],[0,4],[1,3],[1,5]] },
    "C2H2Cl2|1,1-二氯乙烯": { pg: "C2v", mp: "-122.6", bp: "31.6", atoms: [{elem:"C",x:-32,y:0,z:0},{elem:"C",x:32,y:0,z:0},{elem:"Cl",x:-69,y:65,z:0},{elem:"Cl",x:-69,y:-65,z:0},{elem:"H",x:57,y:43,z:0},{elem:"H",x:57,y:-43,z:0}], bonds: [[0,1,"double"],[0,2],[0,3],[1,4],[1,5]] }
});

// C5H12 (戊烷)
addMol("C5H12", "C", "sp³", ["鏈狀/四面體","Chain/Tetra"], "109.5°", "-129.8", "36.1", [], [], {
    "C5H12|正戊烷": { pg: "C2h", mp: "-129.8", bp: "36.1", atoms: [{elem:"C",x:-140,y:-20,z:0},{elem:"C",x:-70,y:20,z:0},{elem:"C",x:0,y:-20,z:0},{elem:"C",x:70,y:20,z:0},{elem:"C",x:140,y:-20,z:0},{elem:"H",x:-140,y:-70,z:0},{elem:"H",x:-175,y:-5,z:25},{elem:"H",x:-175,y:-5,z:-25},{elem:"H",x:-70,y:70,z:0},{elem:"H",x:-70,y:20,z:50},{elem:"H",x:0,y:-70,z:0},{elem:"H",x:0,y:-20,z:-50},{elem:"H",x:70,y:70,z:0},{elem:"H",x:70,y:20,z:50},{elem:"H",x:140,y:-70,z:0},{elem:"H",x:175,y:-5,z:25},{elem:"H",x:175,y:-5,z:-25}], bonds: [[0,1],[1,2],[2,3],[3,4],[0,5],[0,6],[0,7],[1,8],[1,9],[2,10],[2,11],[3,12],[3,13],[4,14],[4,15],[4,16]] },
    "C5H12|異戊烷": { pg: "Cs", mp: "-159.9", bp: "27.8", atoms: [{elem:"C",x:-80,y:-20,z:0},{elem:"C",x:-10,y:20,z:0},{elem:"C",x:60,y:-20,z:0},{elem:"C",x:130,y:20,z:0},{elem:"C",x:-10,y:90,z:0},{elem:"H",x:-80,y:-70,z:0},{elem:"H",x:-115,y:-5,z:25},{elem:"H",x:-115,y:-5,z:-25},{elem:"H",x:-10,y:20,z:50},{elem:"H",x:60,y:-70,z:0},{elem:"H",x:60,y:-20,z:-50},{elem:"H",x:130,y:70,z:0},{elem:"H",x:165,y:-5,z:25},{elem:"H",x:165,y:-5,z:-25},{elem:"H",x:-10,y:140,z:0},{elem:"H",x:-45,y:105,z:25},{elem:"H",x:25,y:105,z:25}], bonds: [[0,1],[1,2],[2,3],[1,4],[0,5],[0,6],[0,7],[1,8],[2,9],[2,10],[3,11],[3,12],[3,13],[4,14],[4,15],[4,16]] },
    "C5H12|新戊烷": { pg: "Td", mp: "-16.5", bp: "9.5", atoms: [{elem:"C",x:0,y:0,z:0},{elem:"C",x:0,y:70,z:0},{elem:"C",x:66,y:-23,z:0},{elem:"C",x:-33,y:-23,z:57},{elem:"C",x:-33,y:-23,z:-57},{elem:"H",x:0,y:120,z:0},{elem:"H",x:-35,y:85,z:35},{elem:"H",x:35,y:85,z:35},{elem:"H",x:116,y:-23,z:0},{elem:"H",x:80,y:-60,z:35},{elem:"H",x:80,y:-60,z:-35},{elem:"H",x:-83,y:-23,z:57},{elem:"H",x:-20,y:27,z:80},{elem:"H",x:-20,y:-73,z:80},{elem:"H",x:-83,y:-23,z:-57},{elem:"H",x:-20,y:27,z:-80},{elem:"H",x:-20,y:-73,z:-80}], bonds: [[0,1],[0,2],[0,3],[0,4],[1,5],[1,6],[1,7],[2,8],[2,9],[2,10],[3,11],[3,12],[3,13],[4,14],[4,15],[4,16]] }
});

// C4H8 (丁烯/環丁烷)
addMol("C4H8", "C", "sp3", ["形狀","Shape"], "N/A", "-185.3", "-6.3", [], [], {
    "C4H8|1-丁烯": { pg: "Cs", mp: "-185.3", bp: "-6.3", atoms: [{elem:"C",x:30,y:-35,z:-15},{elem:"C",x:80,y:20,z:5},{elem:"C",x:-40,y:-20,z:20},{elem:"C",x:-90,y:20,z:-10},{elem:"H",x:25,y:-40,z:-65},{elem:"H",x:50,y:-80,z:5},{elem:"H",x:120,y:10,z:-20},{elem:"H",x:70,y:70,z:-15},{elem:"H",x:90,y:30,z:55},{elem:"H",x:-50,y:-40,z:65},{elem:"H",x:-130,y:25,z:15},{elem:"H",x:-80,y:40,z:-60}], bonds: [[0,1],[0,2],[0,4],[0,5],[1,6],[1,7],[1,8],[2,3,"double"],[2,9],[3,10],[3,11]] },
    "C4H8|順-2-丁烯": { pg: "C2v", mp: "-138.9", bp: "3.7", atoms: [{elem:"C",x:-32,y:-40,z:0},{elem:"C",x:32,y:-40,z:0},{elem:"C",x:-75,y:25,z:0},{elem:"C",x:75,y:25,z:0},{elem:"H",x:-57,y:-85,z:0},{elem:"H",x:57,y:-85,z:0},{elem:"H",x:-115,y:15,z:-30},{elem:"H",x:-55,y:65,z:-20},{elem:"H",x:-90,y:35,z:50},{elem:"H",x:115,y:15,z:30},{elem:"H",x:55,y:65,z:20},{elem:"H",x:90,y:35,z:-50}], bonds: [[0,1,"double"],[0,2],[0,4],[1,3],[1,5],[2,6],[2,7],[2,8],[3,9],[3,10],[3,11]] },
    "C4H8|反-2-丁烯": { pg: "C2h", mp: "-105.5", bp: "0.9", atoms: [{elem:"C",x:-30,y:20,z:0},{elem:"C",x:30,y:-20,z:0},{elem:"C",x:-90,y:-10,z:0},{elem:"C",x:90,y:10,z:0},{elem:"H",x:-20,y:70,z:0},{elem:"H",x:20,y:-70,z:0},{elem:"H",x:-115,y:10,z:40},{elem:"H",x:-95,y:-60,z:0},{elem:"H",x:-115,y:10,z:-40},{elem:"H",x:115,y:-10,z:40},{elem:"H",x:95,y:60,z:0},{elem:"H",x:115,y:-10,z:-40}], bonds: [[0,1,"double"],[0,2],[0,4],[1,3],[1,5],[2,6],[2,7],[2,8],[3,9],[3,10],[3,11]] },
    "C4H8|2-甲基丙烯": { pg: "C2v", mp: "-140.3", bp: "-6.9", atoms: [{elem:"C",x:0,y:10,z:0},{elem:"C",x:-60,y:-25,z:0},{elem:"C",x:60,y:-25,z:0},{elem:"C",x:0,y:73,z:0},{elem:"H",x:-65,y:-60,z:40},{elem:"H",x:-65,y:-60,z:-40},{elem:"H",x:-105,y:5,z:0},{elem:"H",x:65,y:-60,z:-40},{elem:"H",x:65,y:-60,z:40},{elem:"H",x:105,y:5,z:0},{elem:"H",x:43,y:98,z:0},{elem:"H",x:-43,y:98,z:0}], bonds: [[0,1],[0,2],[0,3,"double"],[1,4],[1,5],[1,6],[2,7],[2,8],[2,9],[3,10],[3,11]] },
    "C4H8|環丁烷": { pg: "D2d", mp: "-91", bp: "12.5", atoms: [{elem:"C",x:35,y:35,z:10},{elem:"C",x:-35,y:35,z:-10},{elem:"C",x:-35,y:-35,z:10},{elem:"C",x:35,y:-35,z:-10},{elem:"H",x:45,y:45,z:60},{elem:"H",x:60,y:60,z:-30},{elem:"H",x:-45,y:45,z:-60},{elem:"H",x:-60,y:60,z:30},{elem:"H",x:-60,y:-60,z:30},{elem:"H",x:-45,y:-45,z:60},{elem:"H",x:60,y:-60,z:-30},{elem:"H",x:45,y:-45,z:-60}], bonds: [[0,1],[1,2],[2,3],[3,0],[0,4],[0,5],[1,6],[1,7],[2,8],[2,9],[3,10],[3,11]] },
    "C4H8|甲基環丙烷": { pg: "Cs", mp: "-117.2", bp: "0.7", atoms: [{elem:"C",x:-10,y:0,z:35},{elem:"C",x:50,y:-35,z:0},{elem:"C",x:50,y:35,z:0},{elem:"C",x:-70,y:0,z:-10},{elem:"H",x:-20,y:0,z:85},{elem:"H",x:80,y:-60,z:30},{elem:"H",x:40,y:-60,z:-50},{elem:"H",x:40,y:60,z:-50},{elem:"H",x:80,y:60,z:30},{elem:"H",x:-100,y:-45,z:5},{elem:"H",x:-60,y:0,z:-60},{elem:"H",x:-100,y:45,z:5}], bonds: [[0,1],[0,2],[0,3],[0,4],[1,2],[1,5],[1,6],[2,7],[2,8],[3,9],[3,10],[3,11]] }
});


// C8H10 (二甲苯類) - 已中心化調整
addMol("C8H10", "C", "sp²", ["平面/四面體","Planar/Tetra"], "120°", "-47.8", "139", [], [], {
"C8H10|乙苯|Ethylbenzene": { pg: "Cs", symVectors: { "s": [0,0,1] }, atoms: [{elem:"C",x:0,y:63,z:0},{elem:"C",x:54,y:31,z:0},{elem:"C",x:54,y:-31,z:0},{elem:"C",x:0,y:-63,z:0},{elem:"C",x:-54,y:-31,z:0},{elem:"C",x:-54,y:31,z:0},{elem:"C",x:0,y:129,z:0},{elem:"C",x:66,y:141,z:0},{elem:"H",x:96,y:56,z:0},{elem:"H",x:96,y:-56,z:0},{elem:"H",x:0,y:-112,z:0},{elem:"H",x:-96,y:-56,z:0},{elem:"H",x:-96,y:56,z:0},{elem:"H",x:-39,y:139,z:0},{elem:"H",x:66,y:209,z:0},{elem:"H",x:114,y:101,z:0}], bonds: [[0,1,"double"],[1,2],[2,3,"double"],[3,4],[4,5,"double"],[5,0],[0,6],[1,8],[2,9],[3,10],[4,11],[5,12],[6,7],[6,13],[7,14],[7,15]] },
"C8H10|鄰二甲苯|o-Xylene": { pg: "C2v", symVectors: { "C2": [0,1,0], "sv": [0,0,1], "sv2": [1,0,0] }, atoms: [{elem:"C",x:31,y:54,z:0},{elem:"C",x:-31,y:54,z:0},{elem:"C",x:63,y:0,z:0},{elem:"C",x:-63,y:0,z:0},{elem:"C",x:31,y:-54,z:0},{elem:"C",x:-31,y:-54,z:0},{elem:"C",x:66,y:112,z:0},{elem:"C",x:-66,y:112,z:0},{elem:"H",x:112,y:0,z:0},{elem:"H",x:56,y:-97,z:0},{elem:"H",x:-56,y:-97,z:0},{elem:"H",x:-112,y:0,z:0},{elem:"H",x:115,y:104,z:0},{elem:"H",x:56,y:139,z:40},{elem:"H",x:56,y:139,z:-40},{elem:"H",x:-115,y:104,z:0},{elem:"H",x:-56,y:139,z:40},{elem:"H",x:-56,y:139,z:-40}], bonds: [[0,1],[0,2,"double"],[0,6],[1,3,"double"],[1,7],[2,4],[2,8],[3,5],[3,11],[4,5,"double"],[4,9],[5,10],[6,12],[6,13],[6,14],[7,15],[7,16],[7,17]] },
"C8H10|間二甲苯|m-Xylene": { pg: "C2v", symVectors: { "C2": [0,1,0], "sv": [0,0,1], "sv2": [1,0,0] }, atoms: [{elem:"C",x:-54,y:31.5,z:0},{elem:"C",x:54,y:31.5,z:0},{elem:"C",x:0,y:62.5,z:0},{elem:"C",x:-54,y:-31.5,z:0},{elem:"C",x:54,y:-31.5,z:0},{elem:"C",x:0,y:-62.5,z:0},{elem:"C",x:-112,y:64.5,z:0},{elem:"C",x:112,y:64.5,z:0},{elem:"H",x:0,y:111.5,z:0},{elem:"H",x:-96,y:-56.5,z:0},{elem:"H",x:96,y:-56.5,z:0},{elem:"H",x:0,y:-111.5,z:0},{elem:"H",x:-148,y:38.5,z:0},{elem:"H",x:-126,y:74.5,z:46},{elem:"H",x:-108,y:107.5,z:-46},{elem:"H",x:148,y:38.5,z:0},{elem:"H",x:126,y:74.5,z:46},{elem:"H",x:108,y:107.5,z:-46}], bonds: [[0,2,"double"],[0,3],[0,6],[1,2],[1,4,"double"],[1,7],[2,8],[3,5,"double"],[3,9],[4,5],[4,10],[5,11],[6,12],[6,13],[6,14],[7,15],[7,16],[7,17]] },
"C8H10|對二甲苯|p-Xylene": { pg: "D2h", symVectors: { "C2x": [1,0,0], "C2y": [0,1,0], "C2z": [0,0,1], "sh": [0,0,1], "i": [0,0,0] }, atoms: [{elem:"C",x:63,y:0,z:0},{elem:"C",x:-63,y:0,z:0},{elem:"C",x:31,y:-54,z:0},{elem:"C",x:-31,y:-54,z:0},{elem:"C",x:31,y:54,z:0},{elem:"C",x:-31,y:54,z:0},{elem:"C",x:130,y:0,z:0},{elem:"C",x:-130,y:0,z:0},{elem:"H",x:55,y:-97,z:0},{elem:"H",x:-55,y:-97,z:0},{elem:"H",x:55,y:97,z:0},{elem:"H",x:-55,y:97,z:0},{elem:"H",x:148,y:40,z:0},{elem:"H",x:148,y:-40,z:46},{elem:"H",x:147,y:0,z:-46},{elem:"H",x:-148,y:40,z:0},{elem:"H",x:-148,y:-40,z:46},{elem:"H",x:-147,y:0,z:-46}], bonds: [[0,2,"double"],[0,4],[0,6],[1,3,"double"],[1,5],[1,7],[2,3],[2,8],[3,9],[4,5,"double"],[4,10],[5,11],[6,12],[6,13],[6,14],[7,15],[7,16],[7,17]] }
}, null, null, "D2h");


// C6H12O6 (己糖)
addMol("C6H12O6", "C", "sp³", ["鏈狀/環狀","Chain/Ring"], "109.5°", "146", "dec.", [], [], {
    "C6H12O6|半乳糖": { pg: "C1", mp: "167", bp: "dec.", atoms: [{elem:"O",x:35,y:-42,z:-14},{elem:"O",x:14,y:76,z:-59},{elem:"O",x:-95,y:91,z:1},{elem:"O",x:-127,y:-35,z:16},{elem:"O",x:-39,y:-115,z:-30},{elem:"O",x:154,y:-27,z:13},{elem:"C",x:11,y:63,z:4},{elem:"C",x:-54,y:45,z:20},{elem:"C",x:53,y:10,z:18},{elem:"C",x:-71,y:-15,z:-9},{elem:"C",x:-24,y:-63,z:4},{elem:"C",x:118,y:24,z:-1},{elem:"H",x:24,y:104,z:27},{elem:"H",x:-58,y:42,z:69},{elem:"H",x:53,y:0,z:66},{elem:"H",x:-79,y:-9,z:-57},{elem:"H",x:-24,y:-75,z:51},{elem:"H",x:122,y:31,z:-50},{elem:"H",x:136,y:63,z:22},{elem:"H",x:9,y:39,z:-82},{elem:"H",x:-135,y:77,z:12},{elem:"H",x:-137,y:-73,z:-3},{elem:"H",x:-21,y:-149,z:-10},{elem:"H",x:138,y:-61,z:-8}], bonds: [[0,8],[0,10],[1,6],[1,19],[2,7],[2,20],[3,9],[3,21],[4,10],[4,22],[5,11],[5,23],[6,7],[6,8],[6,12],[7,9],[7,13],[8,11],[8,14],[9,10],[9,15],[10,16],[11,17],[11,18]] },
    "C6H12O6|葡萄糖": { pg: "C1", mp: "146", bp: "dec.", atoms: [{elem:"O",x:-26,y:-57,z:14},{elem:"O",x:-36,y:106,z:-13},{elem:"O",x:88,y:88,z:24},{elem:"O",x:132,y:-29,z:-18},{elem:"O",x:54,y:-124,z:14},{elem:"O",x:-148,y:-48,z:-6},{elem:"C",x:-13,y:51,z:13},{elem:"C",x:53,y:44,z:-6},{elem:"C",x:-51,y:-2,z:-8},{elem:"C",x:76,y:-19,z:11},{elem:"C",x:32,y:-68,z:-9},{elem:"C",x:-116,y:2,z:15},{elem:"H",x:-15,y:56,z:62},{elem:"H",x:58,y:52,z:-55},{elem:"H",x:-53,y:-4,z:-58},{elem:"H",x:85,y:-21,z:59},{elem:"H",x:31,y:-73,z:-58},{elem:"H",x:-117,y:2,z:64},{elem:"H",x:-138,y:43,z:-2},{elem:"H",x:-35,y:102,z:-56},{elem:"H",x:71,y:127,z:12},{elem:"H",x:126,y:-27,z:-62},{elem:"H",x:71,y:-117,z:54},{elem:"H",x:-128,y:-84,z:9}], bonds: [[0,8],[0,10],[1,6],[1,19],[2,7],[2,20],[3,9],[3,21],[4,10],[4,22],[5,11],[5,23],[6,7],[6,8],[6,12],[7,9],[7,13],[8,11],[8,14],[9,10],[9,15],[10,16],[11,17],[11,18]] },
    "C6H12O6|果糖": { pg: "C1", mp: "103", bp: "dec.", atoms: [{elem:"O",x:30,y:-58,z:-2},{elem:"O",x:4,y:106,z:-2},{elem:"O",x:60,y:18,z:67},{elem:"O",x:-117,y:57,z:-1},{elem:"O",x:-92,y:-55,z:-46},{elem:"O",x:148,y:-20,z:-14},{elem:"C",x:-9,y:44,z:-13},{elem:"C",x:45,y:5,z:6},{elem:"C",x:-66,y:26,z:21},{elem:"C",x:-77,y:-42,z:15},{elem:"C",x:-21,y:-77,z:32},{elem:"C",x:100,y:18,z:-33},{elem:"H",x:-18,y:39,z:-62},{elem:"H",x:-61,y:38,z:69},{elem:"H",x:-115,y:-56,z:43},{elem:"H",x:-27,y:-125,z:22},{elem:"H",x:-11,y:-73,z:80},{elem:"H",x:91,y:9,z:-81},{elem:"H",x:115,y:64,z:-28},{elem:"H",x:34,y:119,z:-32},{elem:"H",x:102,y:8,z:73},{elem:"H",x:-118,y:53,z:-45},{elem:"H",x:-131,y:-38,z:-54},{elem:"H",x:134,y:-61,z:-18}], bonds: [[0,7],[0,10],[1,6],[1,19],[2,7],[2,20],[3,8],[3,21],[4,9],[4,22],[5,11],[5,23],[6,7],[6,8],[6,12],[7,11],[8,9],[8,13],[9,10],[9,14],[10,15],[10,16],[11,17],[11,18]] }
});

addMol("C4H4O4", "C", "sp²", ["平面","Planar"], "120°", "", "", [], [], {
    "C4H4O4|順丁烯二酸": { pg: "C2v", mp: "131", bp: "135 (dec)", atoms: [{elem:"C",x:-57,y:-5,z:-2},{elem:"C",x:-29,y:49,z:3},{elem:"C",x:-29,y:-65,z:-6},{elem:"C",x:36,y:62,z:6},{elem:"O",x:19,y:-73,z:31},{elem:"O",x:-49,y:-105,z:-38},{elem:"O",x:70,y:28,z:-31},{elem:"O",x:57,y:101,z:38},{elem:"H",x:-106,y:-6,z:-5},{elem:"H",x:-57,y:89,z:7},{elem:"H",x:33,y:-114,z:23},{elem:"H",x:111,y:40,z:-25}], bonds: [[0,1,"double"],[0,2],[0,8],[1,3],[1,9],[2,4],[2,5,"double"],[3,6],[3,7,"double"],[4,10],[6,11]] },
    "C4H4O4|反丁烯二酸": { pg: "C2h", mp: "287", bp: "290 (subl)", atoms: [{elem:"C",x:-15,y:25,z:-8},{elem:"C",x:15,y:-25,z:8},{elem:"C",x:-80,y:31,z:5},{elem:"C",x:80,y:-31,z:-5},{elem:"O",x:-101,y:84,z:-16},{elem:"O",x:-111,y:-5,z:30},{elem:"O",x:101,y:-84,z:16},{elem:"O",x:111,y:6,z:-30},{elem:"H",x:8,y:61,z:-32},{elem:"H",x:-8,y:-61,z:32},{elem:"H",x:-144,y:85,z:-6},{elem:"H",x:144,y:-85,z:6}], bonds: [[0,1,"double"],[0,2],[0,8],[1,3],[1,9],[2,4],[2,5,"double"],[3,6],[3,7,"double"],[4,10],[6,11]] }
});

addMol("C3H6", "C", "sp3", ["四面體","Tetrahedral"], "109.5", "-185.2", "-47.6", [], [], {
    "C3H6|丙烯": { pg: "Cs", mp: "-185.2", bp: "-47.6", atoms: [{elem:"C",x:54,y:7,z:0},{elem:"C",x:-9,y:-28,z:0},{elem:"C",x:-65,y:9,z:0},{elem:"H",x:79,y:-15,z:-40},{elem:"H",x:54,y:60,z:0},{elem:"H",x:79,y:-15,z:40},{elem:"H",x:-9,y:-78,z:0},{elem:"H",x:-105,y:-20,z:0},{elem:"H",x:-65,y:60,z:0}], bonds: [[0,1],[0,3],[0,4],[0,5],[1,2,"double"],[1,6],[2,7],[2,8]] },
    "C3H6|環丙烷": { pg: "D3h", mp: "-127.6", bp: "-32.9", atoms: [{elem:"C",x:0,y:-40,z:0},{elem:"C",x:-35,y:20,z:0},{elem:"C",x:35,y:20,z:0},{elem:"H",x:0,y:-75,z:-40},{elem:"H",x:0,y:-75,z:40},{elem:"H",x:-65,y:35,z:40},{elem:"H",x:-65,y:35,z:-40},{elem:"H",x:65,y:35,z:-40},{elem:"H",x:65,y:35,z:40}], bonds: [[0,1],[0,2],[0,3],[0,4],[1,2],[1,5],[1,6],[2,7],[2,8]] }
});

addMol("C3H8O", "C", "sp3", ["正四面體","Tetrahedral"], "109.5°", "-126", "97.2", [], [], {
    "C3H8O|1-丙醇": { pg: "Cs", mp: "-126", bp: "97.2", atoms: [{elem:"O",x:115,y:10,z:-5},{elem:"C",x:-5,y:25,z:-5},{elem:"C",x:55,y:-25,z:-5},{elem:"C",x:-70,y:-10,z:-5},{elem:"H",x:0,y:55,z:-45},{elem:"H",x:0,y:55,z:40},{elem:"H",x:50,y:-55,z:-45},{elem:"H",x:50,y:-55,z:40},{elem:"H",x:-75,y:-40,z:-45},{elem:"H",x:-110,y:25,z:-5},{elem:"H",x:-75,y:-40,z:40},{elem:"H",x:120,y:40,z:35}], bonds: [[0,2],[0,11],[1,2],[1,3],[1,4],[1,5],[2,6],[2,7],[3,8],[3,9],[3,10]] },
    "C3H8O|2-丙醇": { pg: "Cs", mp: "-89", bp: "82.3", atoms: [{elem:"O",x:0,y:-80,z:-5},{elem:"C",x:0,y:-15,z:25},{elem:"C",x:-60,y:20,z:0},{elem:"C",x:60,y:20,z:0},{elem:"H",x:0,y:-15,z:75},{elem:"H",x:-60,y:70,z:20},{elem:"H",x:-100,y:-10,z:20},{elem:"H",x:-60,y:20,z:-50},{elem:"H",x:60,y:70,z:20},{elem:"H",x:100,y:-10,z:20},{elem:"H",x:60,y:20,z:-50},{elem:"H",x:0,y:-80,z:-55}], bonds: [[0,1],[0,11],[1,2],[1,3],[1,4],[2,5],[2,6],[2,7],[3,8],[3,9],[3,10]] },
    "C3H8O|甲乙醚": { pg: "Cs", mp: "-113", bp: "7.4", atoms: [{elem:"O",x:0,y:20,z:0},{elem:"C",x:-50,y:-15,z:0},{elem:"C",x:50,y:-15,z:0},{elem:"C",x:95,y:20,z:0},{elem:"H",x:-50,y:-45,z:25},{elem:"H",x:-50,y:-45,z:-25},{elem:"H",x:-85,y:5,z:0},{elem:"H",x:50,y:-45,z:25},{elem:"H",x:50,y:-45,z:-25},{elem:"H",x:95,y:50,z:0},{elem:"H",x:130,y:0,z:25},{elem:"H",x:130,y:0,z:-25}], bonds: [[0,1],[0,2],[1,4],[1,5],[1,6],[2,3],[2,7],[2,8],[3,9],[3,10],[3,11]] }
});

addMol("C3H6Cl2", "C", "sp3", ["正四面體","Tetrahedral"], "109.5°", "-78", "87", [], [], {
    "C3H6Cl2|1,1-二氯丙烷": { pg: "Cs", mp: "-78", bp: "87", atoms: [{elem:"Cl",x:85,y:60,z:-15},{elem:"Cl",x:85,y:-60,z:-15},{elem:"C",x:-15,y:0,z:30},{elem:"C",x:55,y:0,z:20},{elem:"C",x:-55,y:0,z:-35},{elem:"H",x:-30,y:-45,z:55},{elem:"H",x:-30,y:45,z:55},{elem:"H",x:75,y:0,z:70},{elem:"H",x:-45,y:45,z:-65},{elem:"H",x:-105,y:0,z:-25},{elem:"H",x:-45,y:-45,z:-65}], bonds: [[0,3],[1,3],[2,3],[2,4],[2,5],[2,6],[3,7],[4,8],[4,9],[4,10]] },
    "C3H6Cl2|1,2-二氯丙烷": { pg: "C1", mp: "-100.4", bp: "96", atoms: [{elem:"Cl",x:-45,y:-80,z:20},{elem:"Cl",x:95,y:-30,z:30},{elem:"C",x:-25,y:-15,z:-15},{elem:"C",x:45,y:-20,z:-30},{elem:"C",x:-35,y:45,z:30},{elem:"H",x:-50,y:-10,z:-60},{elem:"H",x:65,y:25,z:-55},{elem:"H",x:55,y:-60,z:-60},{elem:"H",x:-15,y:40,z:75},{elem:"H",x:-85,y:50,z:35},{elem:"H",x:-20,y:90,z:10}], bonds: [[0,2],[1,3],[2,3],[2,4],[2,5],[3,6],[3,7],[4,8],[4,9],[4,10]] },
    "C3H6Cl2|1,3-二氯丙烷": { pg: "C2v", mp: "-99.5", bp: "120.4", atoms: [{elem:"Cl",x:90,y:50,z:-15},{elem:"Cl",x:-90,y:50,z:15},{elem:"C",x:0,y:-30,z:0},{elem:"C",x:50,y:10,z:35},{elem:"C",x:-50,y:10,z:-35},{elem:"H",x:-25,y:-65,z:35},{elem:"H",x:25,y:-65,z:-35},{elem:"H",x:85,y:-25,z:55},{elem:"H",x:30,y:40,z:70},{elem:"H",x:-30,y:40,z:-70},{elem:"H",x:-85,y:-25,z:-55}], bonds: [[0,3],[1,4],[2,3],[2,4],[2,5],[2,6],[3,7],[3,8],[4,9],[4,10]] },
    "C3H6Cl2|2,2-二氯丙烷": { pg: "C2v", mp: "-33.8", bp: "69.3", atoms: [{elem:"Cl",x:-60,y:-70,z:0},{elem:"Cl",x:60,y:-70,z:0},{elem:"C",x:0,y:-25,z:0},{elem:"C",x:0,y:20,z:60},{elem:"C",x:0,y:20,z:-60},{elem:"H",x:0,y:-10,z:105},{elem:"H",x:-45,y:50,z:60},{elem:"H",x:45,y:50,z:60},{elem:"H",x:45,y:50,z:-60},{elem:"H",x:0,y:-10,z:-105},{elem:"H",x:-45,y:50,z:-60}], bonds: [[0,2],[1,2],[2,3],[2,4],[3,5],[3,6],[3,7],[4,8],[4,9],[4,10]] }
});

addMol("C3H7Cl", "C", "sp3", ["正四面體","Tetrahedral"], "109.5°", "-122.8", "46.6", [], [], {
    "C3H7Cl|1-氯丙烷": { pg: "Cs", mp: "-122.8", bp: "46.6", atoms: [{elem:"Cl",x:105,y:28,z:-5},{elem:"C",x:-15,y:-22,z:-18},{elem:"C",x:55,y:-30,z:12},{elem:"C",x:-50,y:40,z:5},{elem:"H",x:-12,y:-20,z:-68},{elem:"H",x:-45,y:-65,z:-5},{elem:"H",x:52,y:-32,z:62},{elem:"H",x:78,y:-75,z:-5},{elem:"H",x:-52,y:42,z:55},{elem:"H",x:-100,y:40,z:-15},{elem:"H",x:-25,y:85,z:-12}], bonds: [[0,2],[1,2],[1,3],[1,4],[1,5],[2,6],[2,7],[3,8],[3,9],[3,10]] },
    "C3H7Cl|2-氯丙烷": { pg: "Cs", mp: "-117.2", bp: "35.7", atoms: [{elem:"Cl",x:90,y:0,z:-5},{elem:"C",x:15,y:0,z:20},{elem:"C",x:-20,y:60,z:-5},{elem:"C",x:-20,y:-60,z:-5},{elem:"H",x:14,y:0,z:70},{elem:"H",x:5,y:105,z:15},{elem:"H",x:-65,y:60,z:10},{elem:"H",x:-20,y:60,z:-55},{elem:"H",x:5,y:-105,z:15},{elem:"H",x:-65,y:-60,z:10},{elem:"H",x:-20,y:-60,z:-55}], bonds: [[0,1],[1,2],[1,3],[1,4],[2,5],[2,6],[2,7],[3,8],[3,9],[3,10]] }
});


addMol("C2H4Cl2|二氯乙烷", "C", "sp3", ["正四面體","Tetrahedral"], "109.5°", "-96.9", "57.3", [], [], 
    {
        "C2H4Cl2|1,1-二氯乙烷": { pg: "Cs", mp: "-96.9", bp: "57.3", atoms: [{elem:"Cl",x:66,y:68,z:-9},{elem:"Cl",x:-66,y:68,z:-9},{elem:"C",x:0,y:28,z:15},{elem:"C",x:0,y:-35,z:-10},{elem:"H",x:0,y:27,z:64},{elem:"H",x:0,y:-35,z:-60},{elem:"H",x:-40,y:60,z:5},{elem:"H",x:40,y:60,z:5}], bonds: [[0,2],[1,2],[2,3],[2,4],[3,5],[3,6],[3,7]] },
        "C2H4Cl2|1,2-二氯乙烷": { pg: "C2h", mp: "-35.7", bp: "83.5", atoms: [{elem:"Cl",x:74,y:58,z:6},{elem:"Cl",x:-74,y:58,z:-6},{elem:"C",x:31,y:-6,z:-15},{elem:"C",x:-31,y:-6,z:15},{elem:"H",x:27,y:-6,z:-64},{elem:"H",x:56,y:-46,z:-2},{elem:"H",x:-56,y:-46,z:2},{elem:"H",x:-27,y:-6,z:64}], bonds: [[0,2],[1,3],[2,3],[2,4],[2,5],[3,6],[3,7]] }
    }
);

addMol("C6H4Cl2|二氯苯", "C", "sp²", ["平面 (苯環)","Planar"], "120°", "-17", "180.5", [], [], {
"C6H4Cl2|鄰二氯苯|1,2-Dichlorobenzene": { pg: "C2v", mp: "-17", bp: "180.5", atoms: [{elem:"C",x:-39,y:68,z:0},{elem:"C",x:39,y:68,z:0},{elem:"C",x:78,y:0,z:0},{elem:"C",x:39,y:-68,z:0},{elem:"C",x:-39,y:-68,z:0},{elem:"C",x:-78,y:0,z:0},{elem:"Cl",x:-68,y:117,z:0},{elem:"Cl",x:68,y:117,z:0},{elem:"H",x:135,y:0,z:0},{elem:"H",x:68,y:-117,z:0},{elem:"H",x:-68,y:-117,z:0},{elem:"H",x:-135,y:0,z:0}], bonds: [[0,1,"double"],[1,2],[2,3,"double"],[3,4],[4,5,"double"],[5,0],[0,6],[1,7],[2,8],[3,9],[4,10],[5,11]] },
"C6H4Cl2|間二氯苯|1,3-Dichlorobenzene": { pg: "C2v", mp: "-24.7", bp: "173", atoms: [{elem:"C",x:-68,y:39,z:0},{elem:"C",x:0,y:78,z:0},{elem:"C",x:68,y:39,z:0},{elem:"C",x:68,y:-39,z:0},{elem:"C",x:0,y:-78,z:0},{elem:"C",x:-68,y:-39,z:0},{elem:"Cl",x:-117,y:68,z:0},{elem:"Cl",x:117,y:68,z:0},{elem:"H",x:0,y:135,z:0},{elem:"H",x:117,y:-68,z:0},{elem:"H",x:0,y:-135,z:0},{elem:"H",x:-117,y:-68,z:0}], bonds: [[0,1,"double"],[1,2],[2,3,"double"],[3,4],[4,5,"double"],[5,0],[0,6],[2,7],[1,8],[3,9],[4,10],[5,11]] },
"C6H4Cl2|對二氯苯|1,4-Dichlorobenzene": { pg: "D2h", mp: "53.5", bp: "174", symVectors: { "C2_1": [1,0,0], "C2_2": [0,1,0], "C2_3": [0,0,1], "sh": [0,0,1], "sv1": [1,0,0], "sv2": [0,1,0] }, atoms: [{elem:"C",x:78,y:0,z:0},{elem:"C",x:39,y:-68,z:0},{elem:"C",x:-39,y:-68,z:0},{elem:"C",x:-78,y:0,z:0},{elem:"C",x:-39,y:68,z:0},{elem:"C",x:39,y:68,z:0},{elem:"Cl",x:135,y:0,z:0},{elem:"Cl",x:-135,y:0,z:0},{elem:"H",x:68,y:-117,z:0},{elem:"H",x:-68,y:-117,z:0},{elem:"H",x:-68,y:117,z:0},{elem:"H",x:68,y:117,z:0}], bonds: [[0,1,"double"],[1,2],[2,3,"double"],[3,4],[4,5,"double"],[5,0],[0,6],[3,7],[1,8],[2,9],[4,10],[5,11]] }
});

addMol("C2H2Cl4", "C", "sp3", ["四面體","Tetrahedral"], "109.5°", "", "", [], [], 
{
        "C2H2Cl4|1,1,2,2-四氯乙烷": { pg: "C2h", mp: "-42.3", bp: "146", atoms: [{elem:"Cl",x:-38,y:66,z:59},{elem:"Cl",x:-93,y:-33,z:-9},{elem:"Cl",x:38,y:-66,z:59},{elem:"Cl",x:93,y:33,z:-9},{elem:"C",x:-30,y:17,z:-5},{elem:"C",x:30,y:-17,z:-5},{elem:"H",x:-33,y:46,z:-45},{elem:"H",x:33,y:-46,z:-45}], bonds: [[0,4],[1,4],[2,5],[3,5],[4,5],[4,6],[5,7]] },
        "C2H2Cl4|1,1,1,2-四氯乙烷": { pg: "Cs", mp: "-70.2", bp: "130.2", atoms: [{elem:"Cl",x:27,y:-65,z:-65},{elem:"Cl",x:95,y:-1,z:26},{elem:"Cl",x:29,y:66,z:-63},{elem:"Cl",x:-98,y:0,z:-13},{elem:"C",x:27,y:0,z:-18},{elem:"C",x:-27,y:0,z:25},{elem:"H",x:-27,y:-40,z:54},{elem:"H",x:-26,y:40,z:54}], bonds: [[0,4],[1,4],[2,4],[3,5],[4,5],[5,6],[5,7]] }
    }
);

addMol("C4H6", "C", "sp/sp²/sp³", ["多種異構物", "Isomers"], "Varies", "-108", "10", [], [], 
{
    "C4H6|1,3-丁二烯": { pg: "C2h", mp: "-108.9", bp: "-4.4", atoms: [{elem:"C",x:-27,y:-18,z:0},{elem:"C",x:27,y:18,z:0},{elem:"C",x:-82,y:6,z:0},{elem:"C",x:82,y:-6,z:0},{elem:"H",x:-22,y:-67,z:0},{elem:"H",x:22,y:67,z:0},{elem:"H",x:-122,y:-23,z:0},{elem:"H",x:-90,y:54,z:0},{elem:"H",x:122,y:23,z:0},{elem:"H",x:90,y:-54,z:0}], bonds: [[0,1],[0,2,"double"],[0,4],[1,3,"double"],[1,5],[2,6],[2,7],[3,8],[3,9]] },
    "C4H6|1,2-丁二烯": { pg: "Cs", mp: "-136.2", bp: "10.9", atoms: [{elem:"C",x:-71,y:-13,z:0},{elem:"C",x:-18,y:28,z:0},{elem:"C",x:38,y:11,z:0},{elem:"C",x:94,y:-6,z:0},{elem:"H",x:-98,y:-5,z:-40},{elem:"H",x:-57,y:-60,z:0},{elem:"H",x:-98,y:-5,z:40},{elem:"H",x:-27,y:76,z:0},{elem:"H",x:105,y:-54,z:0},{elem:"H",x:130,y:26,z:0}], bonds: [[0,1],[0,4],[0,5],[0,6],[1,2,"double"],[1,7],[2,3,"double"],[3,8],[3,9]] },
    "C4H6|2-丁炔": { pg: "D3h", mp: "-32.3", bp: "27", atoms: [{elem:"C",x:-93,y:0,z:0},{elem:"C",x:93,y:0,z:0},{elem:"C",x:-27,y:0,z:0},{elem:"C",x:27,y:0,z:0},{elem:"H",x:-110,y:15,z:43},{elem:"H",x:-110,y:-45,z:-9},{elem:"H",x:-110,y:30,z:-35},{elem:"H",x:110,y:46,z:1},{elem:"H",x:110,y:-24,z:39},{elem:"H",x:110,y:-22,z:-40}], bonds: [[0,2],[0,4],[0,5],[0,6],[1,3],[1,7],[1,8],[1,9],[2,3,"triple"]] },
    "C4H6|1-甲基環丙烯": { pg: "Cs", mp: "-", bp: "12", atoms: [{elem:"C",x:53,y:-28,z:0},{elem:"C",x:0,y:8,z:0},{elem:"C",x:51,y:36,z:0},{elem:"C",x:-65,y:3,z:0},{elem:"H",x:65,y:-50,z:42},{elem:"H",x:65,y:-50,z:-42},{elem:"H",x:75,y:77,z:0},{elem:"H",x:-80,y:-22,z:40},{elem:"H",x:-86,y:47,z:0},{elem:"H",x:-80,y:-22,z:-40}], bonds: [[0,1],[0,2],[0,4],[0,5],[1,2,"double"],[1,3],[2,6],[3,7],[3,8],[3,9]] },
    "C4H6|1-丁炔": { pg: "Cs", mp: "-125.7", bp: "8.1", atoms: [{elem:"C",x:10,y:31,z:0},{elem:"C",x:55,y:-20,z:0},{elem:"C",x:-52,y:8,z:0},{elem:"C",x:-103,y:-10,z:0},{elem:"H",x:18,y:59,z:40},{elem:"H",x:18,y:59,z:-40},{elem:"H",x:50,y:-49,z:40},{elem:"H",x:50,y:-49,z:40},{elem:"H",x:101,y:-3,z:0},{elem:"H",x:-148,y:-27,z:0}], bonds: [[0,1],[0,2],[0,4],[0,5],[1,6],[1,7],[1,8],[2,3,"triple"],[3,9]] },
    "C4H6|3-甲基環丙烯": { pg: "Cs", mp: "-", bp: "0", atoms: [{elem:"C",x:5,y:0,z:-24},{elem:"C",x:54,y:-29,z:4},{elem:"C",x:54,y:29,z:4},{elem:"C",x:-55,y:0,z:8},{elem:"H",x:3,y:0,z:-72},{elem:"H",x:76,y:-70,z:16},{elem:"H",x:76,y:70,z:16},{elem:"H",x:-81,y:-40,z:-5},{elem:"H",x:-81,y:40,z:-5},{elem:"H",x:-51,y:0,z:57}], bonds: [[0,1],[0,2],[0,3],[0,4],[1,2,"double"],[1,5],[2,6],[3,7],[3,8],[3,9]] },
    "C4H6|雙環丁烷": { pg: "C2v", mp: "-", bp: "8", atoms: [{elem:"C",x:0,y:-34,z:-12},{elem:"C",x:0,y:34,z:-12},{elem:"C",x:51,y:0,z:15},{elem:"C",x:-51,y:0,z:15},{elem:"H",x:0,y:-53,z:-57},{elem:"H",x:0,y:53,z:-57},{elem:"H",x:93,y:0,z:-9},{elem:"H",x:55,y:0,z:63},{elem:"H",x:-55,y:0,z:63},{elem:"H",x:-93,y:0,z:-9}], bonds: [[0,1],[0,2],[0,3],[0,4],[1,2],[1,3],[1,5],[2,6],[2,7],[3,8],[3,9]] },
    "C4H6|環丁烯": { pg: "C2v", mp: "-135.7", bp: "2.5", atoms: [{elem:"C",x:-26,y:-35,z:0},{elem:"C",x:-26,y:35,z:0},{elem:"C",x:42,y:-30,z:0},{elem:"C",x:42,y:30,z:0},{elem:"H",x:-46,y:-55,z:40},{elem:"H",x:-46,y:-55,z:-40},{elem:"H",x:-46,y:55,z:-40},{elem:"H",x:-46,y:55,z:40},{elem:"H",x:77,y:-64,z:0},{elem:"H",x:77,y:64,z:0}], bonds: [[0,1],[0,2],[0,4],[0,5],[1,3],[1,6],[1,7],[2,3,"double"],[2,8],[3,9]] }
});

addMol("C8H88|Conformer3D_COMPOUND_CID_69667", "C", "sp³", ["幾何形狀","Shape"], "角度", "", "", [{elem:"C",x:39,y:133,z:0},{elem:"C",x:-39,y:133,z:0},{elem:"C",x:34,y:58,z:0},{elem:"C",x:-34,y:58,z:0},{elem:"C",x:71,y:0,z:0},{elem:"C",x:-71,y:0,z:0},{elem:"C",x:35,y:-60,z:0},{elem:"C",x:-35,y:-60,z:0},{elem:"H",x:62,y:156,z:-44},{elem:"H",x:61,y:156,z:44},{elem:"H",x:-62,y:156,z:44},{elem:"H",x:-62,y:156,z:-44},{elem:"H",x:125,y:1,z:0},{elem:"H",x:-125,y:1,z:0},{elem:"H",x:62,y:-107,z:0},{elem:"H",x:-62,y:-107,z:0}], [[0,1],[0,2],[0,8],[0,9],[1,3],[1,10],[1,11],[2,3],[2,4,"double"],[3,5,"double"],[4,6],[4,12],[5,7],[5,13],[6,7,"double"],[6,14],[7,15]], null, null, "C2v");



addMol("C8H8|苯乙烯|Styrene", "C", "sp²", ["平面","Planar"], "120°", "-30.6", "145.2", [{elem:"C",x:18,y:11,z:6},{elem:"C",x:-22,y:59,z:1},{elem:"C",x:-4,y:-48,z:6},{elem:"C",x:-84,y:49,z:-3},{elem:"C",x:-66,y:-58,z:2},{elem:"C",x:-106,y:-10,z:-3},{elem:"C",x:82,y:22,z:10},{elem:"C",x:126,y:-14,z:-11},{elem:"H",x:-6,y:105,z:1},{elem:"H",x:26,y:-86,z:10},{elem:"H",x:-115,y:86,z:-7},{elem:"H",x:-83,y:-104,z:2},{elem:"H",x:-154,y:-18,z:-6},{elem:"H",x:95,y:64,z:31},{elem:"H",x:173,y:-1,z:-5},{elem:"H",x:117,y:-56,z:-34}], [[0,1,"double"],[0,2],[0,6],[1,3],[1,8],[2,4,"double"],[2,9],[3,5,"double"],[3,10],[4,5],[4,11],[5,12],[6,7,"double"],[6,13],[7,14],[7,15]], null, null, "Cs");
addMol("C3H3N|丙烯腈|Acrylonitrile", "C", "sp²/sp", ["平面/直線","Planar/Linear"], "120°/180°", "-83.5", "77.3", [{elem:"N",x:-113,y:10,z:0},{elem:"C",x:-1,y:-22,z:0},{elem:"C",x:43,y:18,z:0},{elem:"C",x:-63,y:-4,z:0},{elem:"H",x:9,y:-70,z:0},{elem:"H",x:90,y:3,z:0},{elem:"H",x:35,y:66,z:0}], [[0,3,"triple"],[1,2,"double"],[1,3],[1,4],[2,5],[2,6]], null, null, "Cs");
addMol("C2H3Cl|氯乙烯", "C", "sp3", ["四面體","Tetrahedral"], "109.5°", "-153.8", "-13.8", [{elem:"Cl",x:-91,y:-12,z:0},{elem:"C",x:-20,y:18,z:0},{elem:"C",x:30,y:-15,z:0},{elem:"H",x:-20,y:67,z:0},{elem:"H",x:74,y:6,z:0},{elem:"H",x:28,y:-64,z:0}], [[0,1],[1,2,"double"],[1,3],[2,4],[2,5]], null, null, "Cs");
addMol("C3H8|丙烷", "C", "sp3", ["正四面體","Tetrahedral"], "109.5°", "-187.7", "-42.1", [{elem:"C",x:0,y:35,z:0},{elem:"C",x:-61,y:-17,z:0},{elem:"C",x:61,y:-17,z:0},{elem:"H",x:0,y:65,z:40},{elem:"H",x:0,y:65,z:-40},{elem:"H",x:-65,y:-55,z:40},{elem:"H",x:-65,y:-55,z:-40},{elem:"H",x:-105,y:15,z:0},{elem:"H",x:105,y:15,z:0},{elem:"H",x:65,y:-55,z:40},{elem:"H",x:65,y:-55,z:-40}], [[0,1],[0,2],[0,3],[0,4],[1,5],[1,6],[1,7],[2,8],[2,9],[2,10]], null, null, "C2v");
addMol("C6H5COOH|苯甲酸|安息香酸", "C", "sp²", ["平面","Planar"], "120°", "122", "249", [{elem:"O",x:118,y:50,z:0},{elem:"O",x:125,y:-52,z:0},{elem:"C",x:28,y:-5,z:0},{elem:"C",x:-5,y:-59,z:0},{elem:"C",x:-2,y:50,z:0},{elem:"C",x:-67,y:-58,z:0},{elem:"C",x:-65,y:51,z:0},{elem:"C",x:-98,y:-3,z:0},{elem:"C",x:93,y:-6,z:0},{elem:"H",x:18,y:-102,z:0},{elem:"H",x:21,y:93,z:0},{elem:"H",x:-93,y:-99,z:0},{elem:"H",x:-89,y:94,z:0},{elem:"H",x:-146,y:-2,z:0},{elem:"H",x:162,y:48,z:0}], [[0,8],[0,14],[1,8,"double"],[2,3,"double"],[2,4],[2,8],[3,5],[3,9],[4,6,"double"],[4,10],[5,7,"double"],[5,11],[6,7],[6,12],[7,13]], null, null, "Cs");
addMol("C7H6O3|鄰羥基苯甲酸|Salicylic Acid|水楊酸|柳酸", "C", "sp²", ["平面","Planar"], "120°", "158.6", "211", [{elem:"O",x:27,y:102,z:-7},{elem:"O",x:116,y:-59,z:-20},{elem:"O",x:126,y:34,z:23},{elem:"C",x:28,y:-5,z:3},{elem:"C",x:-3,y:49,z:-2},{elem:"C",x:-4,y:-59,z:8},{elem:"C",x:-66,y:49,z:-3},{elem:"C",x:-67,y:-59,z:7},{elem:"C",x:-98,y:-5,z:2},{elem:"C",x:93,y:-7,z:3},{elem:"H",x:19,y:-102,z:12},{elem:"H",x:-91,y:91,z:-7},{elem:"H",x:-92,y:-101,z:11},{elem:"H",x:-147,y:-4,z:2},{elem:"H",x:-2,y:135,z:-10},{elem:"H",x:160,y:-59,z:-20}], [[0,4],[0,14],[1,9],[1,15],[2,9,"double"],[3,4],[3,5,"double"],[3,9],[4,6,"double"],[5,7],[5,10],[6,8],[6,11],[7,8,"double"],[7,12],[8,13]], null, null, "Cs");
addMol("C9H8O4|阿斯匹靈|乙醯柳酸|乙醯水楊酸", "C", "sp²", ["平面/苯環","Planar"], "109-120°", "136", "140", [{elem:"O",x:51,y:-18,z:35},{elem:"O",x:-36,y:129,z:-33},{elem:"O",x:31,y:105,z:39},{elem:"O",x:76,y:-29,z:-66},{elem:"C",x:-8,y:-20,z:20},{elem:"C",x:-40,y:32,z:6},{elem:"C",x:-37,y:-76,z:19},{elem:"C",x:-101,y:28,z:-9},{elem:"C",x:-98,y:-80,z:4},{elem:"C",x:-130,y:-27,z:-11},{elem:"C",x:-11,y:90,z:7},{elem:"C",x:91,y:-23,z:-14},{elem:"C",x:154,y:-20,z:8},{elem:"H",x:-13,y:-117,z:30},{elem:"H",x:-127,y:68,z:-20},{elem:"H",x:-121,y:-123,z:3},{elem:"H",x:-177,y:-30,z:-23},{elem:"H",x:163,y:-57,z:39},{elem:"H",x:185,y:-24,z:-31},{elem:"H",x:163,y:23,z:29},{elem:"H",x:-16,y:169,z:-33}], [[0,4],[0,11],[1,10],[1,20],[2,10,"double"],[3,11,"double"],[4,5],[4,6,"double"],[5,7,"double"],[5,10],[6,8],[6,13],[7,9],[7,14],[8,9,"double"],[8,15],[9,16],[11,12],[12,17],[12,18],[12,19]], null, null, "C1");
addMol("C3H4|丙二烯|Allene", "C", "sp (中) / sp² (端)", ["直線軸 (兩端垂直)","Linear Axis (Perpendicular Ends)"], "180° (軸) / 90° (面)", "-136", "-34", [{elem:"C",x:53,y:25,z:-3},{elem:"C",x:0,y:0,z:0},{elem:"C",x:-53,y:-25,z:3},{elem:"H",x:63,y:62,z:-33},{elem:"H",x:90,y:9,z:26},{elem:"H",x:-88,y:-8,z:33},{elem:"H",x:-65,y:-63,z:-26}], [[0,1,"double"],[0,3],[0,4],[1,2,"double"],[2,5],[2,6]], null, null, "D2d");
addMol("H2C2O4|乙二酸|草酸", "C", "sp²", ["平面","Planar"], "120°", "189.5", "365", [{elem:"O",x:57,y:54,z:0},{elem:"O",x:-57,y:-54,z:0},{elem:"O",x:61,y:-49,z:0},{elem:"O",x:-61,y:49,z:0},{elem:"C",x:34,y:-2,z:0},{elem:"C",x:-34,y:2,z:0},{elem:"H",x:101,y:53,z:0},{elem:"H",x:-101,y:-53,z:0}], [[0,4],[0,6],[1,5],[1,7],[2,4,"double"],[3,5,"double"],[4,5]], null, null, "C2h");
addMol("C2H4(OH)2|乙二醇|1,2-乙二醇", "C", "sp³", ["扭轉型/四面體","Gauche/Tetrahedral"], "109.5°", "-12.9", "197.3", [{elem:"O",x:-60,y:36,z:-1},{elem:"O",x:65,y:36,z:12},{elem:"C",x:-29,y:-16,z:19},{elem:"C",x:34,y:-16,z:-8},{elem:"H",x:-54,y:-56,z:5},{elem:"H",x:-27,y:-15,z:68},{elem:"H",x:59,y:-56,z:7},{elem:"H",x:33,y:-15,z:-57},{elem:"H",x:-63,y:34,z:-45},{elem:"H",x:42,y:71,z:-1}], [[0,2],[0,8],[1,3],[1,9],[2,3],[2,4],[2,5],[3,6],[3,7]], null, null, "C2");
addMol("C3H5(OH)3|丙三醇|甘油", "C", "sp³", ["鏈狀/四面體","Chain/Tetrahedral"], "109.5°", "17.8", "290", [{elem:"O",x:4,y:70,z:3},{elem:"O",x:111,y:10,z:-19},{elem:"O",x:-106,y:10,z:-19},{elem:"C",x:3,y:10,z:-20},{elem:"C",x:59,y:-22,z:3},{elem:"C",x:-54,y:-21,z:3},{elem:"H",x:3,y:14,z:-69},{elem:"H",x:62,y:-23,z:52},{elem:"H",x:61,y:-68,z:-14},{elem:"H",x:-56,y:-68,z:-13},{elem:"H",x:-56,y:-20,z:52},{elem:"H",x:4,y:68,z:47},{elem:"H",x:108,y:51,z:-4},{elem:"H",x:-141,y:-11,z:-3}], [[0,3],[0,11],[1,4],[1,12],[2,5],[2,13],[3,4],[3,5],[3,6],[4,7],[4,8],[5,9],[5,10]], null, null, "C1");
addMol("C8H9NO2|乙醯胺基苯酚|普拿疼", "C", "sp²", ["平面/四面體","Planar/Tetra"], "120°", "169", ">250", [{elem:"O",x:187,y:-27,z:0},{elem:"O",x:-103,y:-70,z:0},{elem:"N",x:-57,y:26,z:0},{elem:"C",x:4,y:12,z:0},{elem:"C",x:24,y:-47,z:0},{elem:"C",x:46,y:59,z:0},{elem:"C",x:85,y:-60,z:0},{elem:"C",x:108,y:46,z:0},{elem:"C",x:127,y:-14,z:0},{elem:"C",x:-106,y:-15,z:0},{elem:"C",x:-165,y:18,z:0},{elem:"H",x:-6,y:-86,z:0},{elem:"H",x:32,y:106,z:0},{elem:"H",x:-67,y:70,z:0},{elem:"H",x:99,y:-107,z:0},{elem:"H",x:140,y:82,z:0},{elem:"H",x:-169,y:45,z:-41},{elem:"H",x:-202,y:-15,z:1},{elem:"H",x:-168,y:47,z:40},{elem:"H",x:192,y:-70,z:0}], [[0,8],[0,19],[1,9,"double"],[2,3],[2,9],[2,13],[3,4,"double"],[3,5],[4,6],[4,11],[5,7,"double"],[5,12],[6,8,"double"],[6,14],[7,8],[7,15],[9,10],[10,16],[10,17],[10,18]], null, null, "Cs");
addMol("C7H8|甲苯|Toluene", "C", "sp²", ["平面/四面體","Mixed"], "120°", "-95", "110.6", [{elem:"C",x:-62.5,y:0,z:0},{elem:"C",x:-31.5,y:-54,z:0},{elem:"C",x:-31.5,y:54,z:0},{elem:"C",x:-129.5,y:0,z:0},{elem:"C",x:31.5,y:-54,z:0},{elem:"C",x:31.5,y:54,z:0},{elem:"C",x:62.5,y:0,z:0},{elem:"H",x:-55.5,y:-97,z:0},{elem:"H",x:-55.5,y:97,z:0},{elem:"H",x:-147.5,y:40,z:23},{elem:"H",x:-147.5,y:-40,z:23},{elem:"H",x:-146.5,y:0,z:-46},{elem:"H",x:55.5,y:-97,z:0},{elem:"H",x:55.5,y:97,z:0},{elem:"H",x:111.5,y:0,z:0}], [[0,1,"double"],[0,2],[0,3],[1,4],[1,7],[2,5,"double"],[2,8],[3,9],[3,10],[3,11],[4,6,"double"],[4,12],[5,6],[5,13],[6,14]], null, null, "Cs");addMol("C10H8|萘|Naphthalene", "C", "sp²", ["平面","Planar"], "120°", "80.2", "218", [{elem:"C",x:0,y:-32,z:0},{elem:"C",x:0,y:32,z:0},{elem:"C",x:55,y:-63,z:0},{elem:"C",x:55,y:63,z:0},{elem:"C",x:-55,y:-63,z:0},{elem:"C",x:-55,y:63,z:0},{elem:"C",x:109,y:-31,z:0},{elem:"C",x:109,y:31,z:0},{elem:"C",x:-109,y:-31,z:0},{elem:"C",x:-109,y:31,z:0},{elem:"H",x:56,y:-112,z:0},{elem:"H",x:56,y:112,z:0},{elem:"H",x:-56,y:-112,z:0},{elem:"H",x:-56,y:112,z:0},{elem:"H",x:152,y:-56,z:0},{elem:"H",x:152,y:56,z:0},{elem:"H",x:-152,y:-56,z:0},{elem:"H",x:-152,y:56,z:0}], [[0,1],[0,2,"double"],[0,4],[1,3,"double"],[1,5],[2,6],[2,10],[3,7],[3,11],[4,8,"double"],[4,12],[5,9,"double"],[5,13],[6,7,"double"],[6,14],[7,15],[8,9],[8,16],[9,17]], null, null, "D2h");
addMol("C6H6|苯|Benzene", "C", "sp²", ["平面","Planar"], "120°", "5.5", "80.1", [{elem:"C",x:63,y:0,z:0},{elem:"C",x:31.5,y:54.6,z:0},{elem:"C",x:-31.5,y:54.6,z:0},{elem:"C",x:-63,y:0,z:0},{elem:"C",x:-31.5,y:-54.6,z:0},{elem:"C",x:31.5,y:-54.6,z:0},{elem:"H",x:112,y:0,z:0},{elem:"H",x:56,y:97,z:0},{elem:"H",x:-56,y:97,z:0},{elem:"H",x:-112,y:0,z:0},{elem:"H",x:-56,y:-97,z:0},{elem:"H",x:56,y:-97,z:0}], [[0,1,"double"],[1,2],[2,3,"double"],[3,4],[4,5,"double"],[5,0],[0,6],[1,7],[2,8],[3,9],[4,10],[5,11]], null, null, "D6h");
addMol("CH3COOCH3|乙酸甲酯|Methyl Acetate", "C", "sp²", ["平面/四面體","Mixed"], "120°", "-98", "56.9", [{elem:"C",x:-35,y:0,z:0},{elem:"C",x:10,y:0,z:0},{elem:"O",x:10,y:50,z:0},{elem:"O",x:50,y:-30,z:0},{elem:"C",x:90,y:-30,z:0},{elem:"H",x:-35,y:-45,z:0},{elem:"H",x:-65,y:25,z:25},{elem:"H",x:-65,y:25,z:-25},{elem:"H",x:90,y:-70,z:0},{elem:"H",x:120,y:-5,z:25},{elem:"H",x:120,y:-5,z:-25}], [[0,1],[1,2,"double"],[1,3],[3,4],[0,5],[0,6],[0,7],[4,8],[4,9],[4,10]], null, null, "Cs");
addMol("CH3COOC2H5|乙酸乙酯|Ethyl Acetate", "C", "sp²", ["平面/四面體","Mixed"], "120°", "-83.6", "77.1", [{elem:"C",x:-35,y:0,z:0},{elem:"C",x:10,y:0,z:0},{elem:"O",x:10,y:50,z:0},{elem:"O",x:50,y:-30,z:0},{elem:"C",x:90,y:-30,z:0},{elem:"C",x:125,y:15,z:0},{elem:"H",x:-35,y:-45,z:0},{elem:"H",x:-65,y:25,z:25},{elem:"H",x:-65,y:25,z:-25},{elem:"H",x:90,y:-70,z:0},{elem:"H",x:115,y:-50,z:30},{elem:"H",x:125,y:55,z:0},{elem:"H",x:155,y:-5,z:25},{elem:"H",x:155,y:-5,z:-25}], [[0,1],[1,2,"double"],[1,3],[3,4],[4,5],[0,6],[0,7],[0,8],[4,9],[4,10],[5,11],[5,12],[5,13]], null, null, "Cs");
addMol("CH3NHCH3|二甲胺|Dimethylamine", "N", "sp³", ["角錐","Pyramidal"], "107°", "-92.2", "6.9", [{elem:"N",x:0,y:0,z:0},{elem:"H",x:0,y:30,z:0},{elem:"C",x:40,y:-25,z:25},{elem:"C",x:-40,y:-25,z:25},{elem:"H",x:40,y:-25,z:65},{elem:"H",x:70,y:-50,z:5},{elem:"H",x:60,y:10,z:15},{elem:"H",x:-40,y:-25,z:65},{elem:"H",x:-70,y:-50,z:5},{elem:"H",x:-60,y:10,z:15}], [[0,1],[0,2],[0,3],[2,4],[2,5],[2,6],[3,7],[3,8],[3,9]], null, null, "Cs");
addMol("N(CH3)3|三甲胺|Trimethylamine", "N", "sp³", ["角錐","Pyramidal"], "108°", "-117.2", "2.9", [{elem:"N",x:0,y:15,z:0},{elem:"C",x:0,y:-25,z:35},{elem:"C",x:35,y:-25,z:-25},{elem:"C",x:-35,y:-25,z:-25},{elem:"H",x:0,y:-25,z:75},{elem:"H",x:30,y:-55,z:35},{elem:"H",x:-30,y:-55,z:35},{elem:"H",x:35,y:-25,z:-65},{elem:"H",x:65,y:-55,z:-25},{elem:"H",x:65,y:5,z:-25},{elem:"H",x:-35,y:-25,z:-65},{elem:"H",x:-65,y:-55,z:-25},{elem:"H",x:-65,y:5,z:-25}], [[0,1],[0,2],[0,3],[1,4],[1,5],[1,6],[2,7],[2,8],[2,9],[3,10],[3,11],[3,12]], null, null, "C3v");
addMol("HCONH2|甲醯胺|Formamide", "C", "sp²", ["平面","Planar"], "120°", "2.6", "210", [{elem:"C",x:0,y:0,z:0},{elem:"O",x:0,y:50,z:0},{elem:"N",x:45,y:-30,z:0},{elem:"H",x:-40,y:-25,z:0},{elem:"H",x:45,y:-65,z:0},{elem:"H",x:80,y:-10,z:0}], [[0,1,"double"],[0,2],[0,3],[2,4],[2,5]], null, null, "Cs");
addMol("CH3CONH2|乙醯胺|Acetamide", "C", "sp²", ["平面","Planar"], "120°", "82.3", "221.2", [{elem:"C",x:-40,y:0,z:0},{elem:"C",x:5,y:0,z:0},{elem:"O",x:5,y:50,z:0},{elem:"N",x:50,y:-30,z:0},{elem:"H",x:-40,y:-45,z:0},{elem:"H",x:-70,y:25,z:25},{elem:"H",x:-70,y:25,z:-25},{elem:"H",x:50,y:-65,z:0},{elem:"H",x:85,y:-10,z:0}], [[0,1],[1,2,"double"],[1,3],[0,4],[0,5],[0,6],[3,7],[3,8]], null, null, "Cs");
addMol("P4|白磷|黃磷|White Phosphorus", "P", "sp³", ["正四面體 (籠狀)","Tetrahedral Cage"], "60°", "44.1", "280.5", [{elem:"P",x:50,y:50,z:50}, {elem:"P",x:50,y:-50,z:-50}, {elem:"P",x:-50,y:50,z:-50}, {elem:"P",x:-50,y:-50,z:50}], [[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]], null, null, "Td");
addMol("P4O6|六氧化四磷|Phosphorus Hexoxide", "P", "sp³", ["籠狀 (類金剛烷)","Cage"], "100°", "23.8", "173.1", [{elem:"P",x:55,y:55,z:55}, {elem:"P",x:55,y:-55,z:-55}, {elem:"P",x:-55,y:55,z:-55}, {elem:"P",x:-55,y:-55,z:55}, {elem:"O",x:85,y:0,z:0}, {elem:"O",x:-85,y:0,z:0}, {elem:"O",x:0,y:85,z:0}, {elem:"O",x:0,y:-85,z:0}, {elem:"O",x:0,y:0,z:85}, {elem:"O",x:0,y:0,z:-85}], [[0,4],[0,6],[0,8], [1,4],[1,7],[1,9], [2,5],[2,6],[2,9], [3,5],[3,7],[3,8]], null, null, "Td");
addMol("P4O10|十氧化四磷|Phosphorus Pentoxide", "P", "sp³", ["籠狀","Cage"], "102°/123°", "340 (昇華)", "-", [{elem:"P",x:55,y:55,z:55}, {elem:"P",x:55,y:-55,z:-55}, {elem:"P",x:-55,y:55,z:-55}, {elem:"P",x:-55,y:-55,z:55}, {elem:"O",x:85,y:0,z:0}, {elem:"O",x:-85,y:0,z:0}, {elem:"O",x:0,y:85,z:0}, {elem:"O",x:0,y:-85,z:0}, {elem:"O",x:0,y:0,z:85}, {elem:"O",x:0,y:0,z:-85}, {elem:"O",x:95,y:95,z:95}, {elem:"O",x:95,y:-95,z:-95}, {elem:"O",x:-95,y:95,z:-95}, {elem:"O",x:-95,y:-95,z:95}], [[0,4],[0,6],[0,8], [1,4],[1,7],[1,9], [2,5],[2,6],[2,9], [3,5],[3,7],[3,8], [0,10,"double"], [1,11,"double"], [2,12,"double"], [3,13,"double"]], null, null, "Td");
addMol("S8|斜方硫|單斜硫|硫磺", "S", "sp3", ["皇冠型 (環狀)","Crown Ring"], "108°", "115.2", "444.6", [{elem:"S",x:0,y:-105,z:-25},{elem:"S",x:-75,y:-75,z:25},{elem:"S",x:75,y:-75,z:25},{elem:"S",x:-105,y:0,z:-25},{elem:"S",x:105,y:0,z:-25},{elem:"S",x:-75,y:75,z:25},{elem:"S",x:75,y:75,z:25},{elem:"S",x:0,y:105,z:-25}], [[0,1],[0,2],[1,3],[2,4],[3,5],[4,6],[5,7],[6,7]], null, null, "D4d");
addMol("H2O2|過氧化氫|雙氧水", "O", "sp³", ["書本型","Open Book"], "111°", "-0.4", "150.2", [{elem:"O",x:-17,y:5,z:-17},{elem:"O",x:17,y:-5,z:17},{elem:"H",x:-55,y:35,z:-22},{elem:"H",x:55,y:-35,z:22}], [[0,1],[0,2],[1,3]], null, null, "C2");
addMol("H2S2|二硫化氫", "S", "sp³", ["書本型","Open Book"], "92°", "-89.6", "70.7", [{elem:"S",x:-20,y:5,z:-20},{elem:"S",x:20,y:-5,z:20},{elem:"H",x:-55,y:35,z:-25},{elem:"H",x:55,y:-35,z:25}], [[0,1],[0,2],[1,3]], null, null, "C2");
addMol("S2Cl2|二氯化二硫|二氯化硫", "S", "sp³", ["書本型","Open Book"], "103°", "-77", "138", [{elem:"S",x:0,y:25,z:0},{elem:"S",x:0,y:-25,z:0},{elem:"Cl",x:65,y:60,z:50},{elem:"Cl",x:-65,y:-60,z:50}], [[0,1],[0,2],[1,3]], null, null, "C2");
addMol("N2O3|三氧化二氮", "N", "sp²+sp²", ["平面","Planar"], "120°", "-100.7", "3.5", [{elem:"N",x:-32,y:0,z:0},{elem:"N",x:32,y:0,z:0},{elem:"O",x:-67,y:50,z:0},{elem:"O",x:67,y:50,z:0},{elem:"O",x:67,y:-50,z:0}], [[0,1],[0,2,"double"],[1,3,"double"],[1,4]], null, null, "Cs");
addMol("N2O4|四氧化二氮", "N", "sp²", ["平面","Planar"], "120°", "-11.2", "21.2", [{elem:"N",x:-35,y:0,z:0,lpCount:0},{elem:"N",x:35,y:0,z:0,lpCount:0},{elem:"O",x:-75,y:55,z:0},{elem:"O",x:-75,y:-55,z:0},{elem:"O",x:75,y:55,z:0},{elem:"O",x:75,y:-55,z:0}], [[0,1],[0,2,"double"],[0,3],[1,4,"double"],[1,5]], null, null, "D2h");
addMol("N2O5|五氧化二氮", "N", "sp²", ["非平面","V-shape"], "120°", "30 (昇華)", "47 (分解)", [{elem:"O",x:0,y:30,z:0},{elem:"N",x:-50,y:-15,z:0},{elem:"N",x:50,y:-15,z:0},{elem:"O",x:-85,y:30,z:0},{elem:"O",x:-85,y:-60,z:0},{elem:"O",x:85,y:30,z:0},{elem:"O",x:85,y:-60,z:0}], [[0,1],[0,2],[1,3,"double"],[1,4],[2,5,"double"],[2,6]], null, null, "C2");
addMol("C2H5Cl|氯乙烷", "C", "sp3", ["四面體","Tetrahedral"], "109.5", "-138.7", "12.3", [{elem:"Cl",x:100,y:-16,z:0},{elem:"C",x:25,y:28,z:0},{elem:"C",x:-35,y:-20,z:0},{elem:"H",x:25,y:56,z:-40},{elem:"H",x:25,y:56,z:40},{elem:"H",x:-78,y:5,z:0},{elem:"H",x:-35,y:-55,z:40},{elem:"H",x:-35,y:-55,z:-40}], [[0,1],[1,2],[1,3],[1,4],[2,5],[2,6],[2,7]], null, null, "Cs");
addMol("CHOCHO|乙二醛|Glyoxal", "C", "sp²", ["平面","Planar"], "120°", "15", "50.4", [{elem:"C",x:-28,y:0,z:0},{elem:"C",x:28,y:0,z:0},{elem:"O",x:-55,y:45,z:0},{elem:"O",x:55,y:-45,z:0},{elem:"H",x:-45,y:-40,z:0},{elem:"H",x:45,y:40,z:0}], [[0,1],[0,2,"double"],[1,3,"double"],[0,4],[1,5]], null, null, "C2h");
addMol("CH3OH|甲醇|木精", "C", "sp³", ["四面體","Tetrahedral"], "109.5°", "-97.6", "64.7", [{elem:"C",x:0,y:0,z:0}, {elem:"O",x:0,y:70,z:0,lpCount:2}, {elem:"H",x:0,y:-50,z:0}, {elem:"H",x:45,y:15,z:30}, {elem:"H",x:-45,y:15,z:30}, {elem:"H",x:50,y:85,z:0}], [[0,1],[0,2],[0,3],[0,4],[1,5]], null, null, "Cs");
addMol("C2H5OH|乙醇|Ethanol|酒精", "C", "sp3", ["四面體","Tetrahedral"], "109.5°", "-114.1", "78.2", [{elem:"O",x:-69,y:-14,z:-4},{elem:"C",x:-19,y:25,z:-4},{elem:"C",x:38,y:-13,z:-4},{elem:"H",x:-21,y:54,z:36},{elem:"H",x:-21,y:53,z:-44},{elem:"H",x:78,y:16,z:-5},{elem:"H",x:39,y:-42,z:-43},{elem:"H",x:40,y:-41,z:36},{elem:"H",x:-67,y:-38,z:32}], [[0,1],[0,8],[1,2],[1,3],[1,4],[2,5],[2,6],[2,7]], null, null, "Cs");
addMol("C10H15N|甲基苯丙胺|Methamphetamine", "C,N", "sp²/sp³", ["苯環平面","側鏈四面體"], "120°/109.5°", "3", "212", [{elem:"N",x:-126,y:-6,z:0,lpCount:1},{elem:"C",x:-63,y:-12,z:-18},{elem:"C",x:-22,y:11,z:33},{elem:"C",x:44,y:8,z:18},{elem:"C",x:-51,y:-79,z:-32},{elem:"C",x:71,y:57,z:-10},{elem:"C",x:76,y:-43,z:33},{elem:"C",x:-140,y:56,z:15},{elem:"C",x:132,y:54,z:-23},{elem:"C",x:137,y:-45,z:19},{elem:"C",x:165,y:3,z:-9},{elem:"H",x:-55,y:14,z:-59},{elem:"H",x:-32,y:58,z:44},{elem:"H",x:-31,y:-14,z:75},{elem:"H",x:-7,y:-85,z:-53},{elem:"H",x:-55,y:-107,z:8},{elem:"H",x:-84,y:-95,z:-65},{elem:"H",x:-153,y:-19,z:-35},{elem:"H",x:46,y:97,z:-21},{elem:"H",x:55,y:-81,z:55},{elem:"H",x:-189,y:62,z:16},{elem:"H",x:-125,y:68,z:60},{elem:"H",x:-123,y:88,z:-18},{elem:"H",x:154,y:92,z:-45},{elem:"H",x:163,y:-85,z:30},{elem:"H",x:213,y:1,z:-20}], [[0,1],[0,7],[0,17],[1,2],[1,4],[1,11],[2,3],[2,12],[2,13],[3,5,"double"],[3,6],[4,14],[4,15],[4,16],[5,8],[5,18],[6,9,"double"],[6,19],[7,20],[7,21],[7,22],[8,10,"double"],[8,23],[9,10],[9,24],[10,25]], null, null, "C1");


addMol("aaa", "C", "sp3", ["形狀","Shape"], "角度", "", "", [{elem:C,x:20,y:12,z:6},{elem:C,x:-24,y:66,z:1},{elem:C,x:-4,y:-53,z:7},{elem:C,x:-93,y:54,z:-3},{elem:C,x:-73,y:-65,z:2},{elem:C,x:-117,y:-11,z:-3},{elem:C,x:91,y:24,z:11},{elem:C,x:140,y:-16,z:-12},{elem:H,x:-6,y:117,z:1},{elem:H,x:29,y:-96,z:11},{elem:H,x:-127,y:96,z:-7},{elem:H,x:-92,y:-116,z:2},{elem:H,x:-171,y:-20,z:-7},{elem:H,x:106,y:71,z:35},{elem:H,x:192,y:-1,z:-6},{elem:H,x:130,y:-62,z:-38}], [[0,1,'double'],[0,2],[0,6],[1,3],[1,8],[2,4,'double'],[2,9],[3,5,'double'],[3,10],[4,5],[4,11],[5,12],[6,7,'double'],[6,13],[7,14],[7,15]]);
addMol("bbb", "C", "sp3", ["形狀","Shape"], "角度", "", "", [{elem:C,x:91,y:39,z:0},{elem:C,x:91,y:-39,z:0},{elem:C,x:16,y:34,z:0},{elem:C,x:16,y:-34,z:0},{elem:C,x:-42,y:71,z:0},{elem:C,x:-42,y:-71,z:0},{elem:C,x:-102,y:35,z:0},{elem:C,x:-102,y:-35,z:0},{elem:H,x:114,y:62,z:44},{elem:H,x:114,y:61,z:-44},{elem:H,x:114,y:-62,z:-44},{elem:H,x:114,y:-62,z:44},{elem:H,x:-42,y:125,z:0},{elem:H,x:-42,y:-125,z:0},{elem:H,x:-149,y:62,z:0},{elem:H,x:-149,y:-62,z:0}], [[0,1],[0,2],[0,8],[0,9],[1,3],[1,10],[1,11],[2,3],[2,4,'double'],[3,5,'double'],[4,6],[4,12],[5,7],[5,13],[6,7,'double'],[6,14],[7,15]]);







// --- 13.簡單離子化合物
addMol("KCl|氯化鉀", "K", "-", "-", "-", "770", "1420", [{elem:"K",x:-50,y:0,z:0,r:22,lpCount:0}, {elem:"Cl",x:50,y:0,z:0,r:35,lpCount:0}], [[0, 1, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">🧂 物質性質</div><div class="info-body"><strong>氯化鉀 (KCl)</strong><br>白色結晶固體，外觀與食鹽相似。它是鉀肥的主要成分，對植物生長至關重要。</div></div>');
addMol("KI|碘化鉀", "K", "-", "-", "-", "681", "1330", [{elem:"K",x:-55,y:0,z:0,r:22,lpCount:0}, {elem:"I",x:55,y:0,z:0,r:40,lpCount:0}], [[0, 1, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">🧂 物質性質</div><div class="info-body"><strong>碘化鉀 (KI)</strong><br>白色晶體，易溶於水。常添加於食鹽中作為碘的來源。</div></div>');
addMol("KBr|溴化鉀", "K", "-", "-", "-", "734", "1435", [{elem:"K",x:-50,y:0,z:0,r:22,lpCount:0}, {elem:"Br",x:50,y:0,z:0,r:38,lpCount:0}], [[0, 1, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">🧂 物質性質</div><div class="info-body"><strong>溴化鉀 (KBr)</strong><br>白色結晶，對紅外光透明，常用來製作光譜分析的樣品鹽片。</div></div>');
addMol("NaF|氟化鈉", "Na", "-", "-", "-", "993", "1704", [{elem:"Na",x:-40,y:0,z:0,r:20,lpCount:0}, {elem:"F",x:40,y:0,z:0,r:25,lpCount:0}], [[0, 1, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">🧂 物質性質</div><div class="info-body"><strong>氟化鈉 (NaF)</strong><br>牙膏中常見的添加劑，能提供氟離子以強化牙齒琺瑯質。</div></div>');
addMol("LiF|氟化鋰", "Li", "-", "-", "-", "845", "1676", [{elem:"Li",x:-40,y:0,z:0,r:15,lpCount:0}, {elem:"F",x:40,y:0,z:0,r:25,lpCount:0}], [[0, 1, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">🧂 物質性質</div><div class="info-body"><strong>氟化鋰 (LiF)</strong><br>在紫外線區域具有極佳的穿透性，常用於光學透鏡材料。</div></div>');
addMol("MgO|氧化鎂|苦土", "Mg", "-", "-", "-", "2852", "3600", [{elem:"Mg",x:-40,y:0,z:0,r:18,lpCount:0}, {elem:"O",x:40,y:0,z:0,r:25,lpCount:0}], [[0, 1, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">🧱 物質性質</div><div class="info-body"><strong>氧化鎂 (MgO)</strong><br>熔點極高，是優良的耐火材料。</div></div>');
addMol("CaO|氧化鈣|生石灰", "Ca", "-", "-", "-", "2572", "2850", [{elem:"Ca",x:-45,y:0,z:0,r:22,lpCount:0}, {elem:"O",x:45,y:0,z:0,r:25,lpCount:0}], [[0, 1, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">🧱 物質性質</div><div class="info-body"><strong>氧化鈣 (CaO)</strong><br>俗稱生石灰，遇水放熱生成熟石灰，是常用的乾燥劑。</div></div>');
addMol("BaO|氧化鋇", "Ba", "-", "-", "-", "1923", "2000", [{elem:"Ba",x:-50,y:0,z:0,r:28,lpCount:0}, {elem:"O",x:50,y:0,z:0,r:25,lpCount:0}], [[0, 1, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">🧱 物質性質</div><div class="info-body"><strong>氧化鋇 (BaO)</strong><br>用於玻璃工業增加折射率。</div></div>');
addMol("ZnO|氧化鋅|鋅白", "Zn", "-", "-", "-", "1975", "-", [{elem:"Zn",x:-40,y:0,z:0,r:18,lpCount:0}, {elem:"O",x:40,y:0,z:0,r:25,lpCount:0}], [[0, 1, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">🎨 物質性質</div><div class="info-body"><strong>氧化鋅 (ZnO)</strong><br>俗稱鋅白，具有紫外線遮蔽能力，用於防曬乳與橡膠工業。</div></div>');
addMol("CuO|氧化銅", "Cu", "-", "-", "-", "1326", "-", [{elem:"Cu",x:-40,y:0,z:0,r:18,lpCount:0}, {elem:"O",x:40,y:0,z:0,r:25,lpCount:0}], [[0, 1, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">⚫ 物質性質</div><div class="info-body"><strong>氧化銅 (CuO)</strong><br>黑色固體，用於製造顏料與有機分析。</div></div>');
addMol("AgCl|氯化銀", "Ag", "-", "-", "-", "455", "1550", [{elem:"Ag",x:-45,y:0,z:0,r:22,lpCount:0}, {elem:"Cl",x:45,y:0,z:0,r:35,lpCount:0}], [[0, 1, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">📷 物質性質</div><div class="info-body"><strong>氯化銀 (AgCl)</strong><br>白色沈澱，見光分解產生黑色的銀，曾用於攝影底片。</div></div>');
addMol("AgBr|溴化銀", "Ag", "-", "-", "-", "432", "1502", [{elem:"Ag",x:-48,y:0,z:0,r:22,lpCount:0}, {elem:"Br",x:48,y:0,z:0,r:38,lpCount:0}], [[0, 1, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">📷 物質性質</div><div class="info-body"><strong>溴化銀 (AgBr)</strong><br>淺黃色固體，感光性強，傳統攝影底片的主要成分。</div></div>');
addMol("AgI|碘化銀", "Ag", "-", "-", "-", "558", "1506", [{elem:"Ag",x:-50,y:0,z:0,r:22,lpCount:0}, {elem:"I",x:50,y:0,z:0,r:40,lpCount:0}], [[0, 1, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">🌧️ 物質性質</div><div class="info-body"><strong>碘化銀 (AgI)</strong><br>黃色固體，晶體結構似冰，用於人造降雨的晶種。</div></div>');
addMol("NaH|氫化鈉", "Na", "-", "-", "-", "800", "分解", [{elem:"Na",x:-40,y:0,z:0,r:20,lpCount:0}, {elem:"H",x:40,y:0,z:0,r:15,lpCount:0}], [[0, 1, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">🧪 物質性質</div><div class="info-body"><strong>氫化鈉 (NaH)</strong><br>含氫負離子 (H⁻) 的強還原劑。</div></div>');
addMol("HgS|硫化汞|硃砂", "Hg", "-", "-", "-", "583", "昇華", [{elem:"Hg",x:-45,y:0,z:0,r:25,lpCount:0}, {elem:"S",x:45,y:0,z:0,r:30,lpCount:0}], [[0, 1, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">🎨 物質性質</div><div class="info-body"><strong>硫化汞 (HgS)</strong><br>天然硃砂，鮮紅色，古代顏料與煉丹原料。</div></div>');
addMol("MgCl2|氯化鎂", "Mg", "-", "-", "-", "714", "1412", [{elem:"Mg",x:0,y:0,z:0,r:18,lpCount:0}, {elem:"Cl",x:-85,y:0,z:0,r:35,lpCount:0}, {elem:"Cl",x:85,y:0,z:0,r:35,lpCount:0}], [[0, 1, "ionic_thin"], [0, 2, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">🧂 物質性質</div><div class="info-body"><strong>氯化鎂 (MgCl₂)</strong><br>苦滷的主要成分，豆腐凝固劑。</div></div>');
addMol("CaCl2|氯化鈣", "Ca", "-", "-", "-", "772", "1935", [{elem:"Ca",x:0,y:0,z:0,r:22,lpCount:0}, {elem:"Cl",x:-90,y:0,z:0,r:35,lpCount:0}, {elem:"Cl",x:90,y:0,z:0,r:35,lpCount:0}], [[0, 1, "ionic_thin"], [0, 2, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">🧂 物質性質</div><div class="info-body"><strong>氯化鈣 (CaCl₂)</strong><br>強吸濕性，常用乾燥劑與融雪劑。</div></div>');
addMol("CaF2|氟化鈣|螢石", "Ca", "-", "-", "-", "1418", "2533", [{elem:"Ca",x:0,y:0,z:0,r:22,lpCount:0}, {elem:"F",x:-80,y:0,z:0,r:25,lpCount:0}, {elem:"F",x:80,y:0,z:0,r:25,lpCount:0}], [[0, 1, "ionic_thin"], [0, 2, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">💎 物質性質</div><div class="info-body"><strong>氟化鈣 (CaF₂)</strong><br>螢石，製備 HF 的原料，也用於光學鏡頭。</div></div>');
addMol("BaCl2|氯化鋇", "Ba", "-", "-", "-", "962", "1560", [{elem:"Ba",x:0,y:0,z:0,r:28,lpCount:0}, {elem:"Cl",x:-95,y:0,z:0,r:35,lpCount:0}, {elem:"Cl",x:95,y:0,z:0,r:35,lpCount:0}], [[0, 1, "ionic_thin"], [0, 2, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">🧂 物質性質</div><div class="info-body"><strong>氯化鋇 (BaCl₂)</strong><br>檢驗硫酸根的試劑，劇毒，燃燒呈黃綠色火焰。</div></div>');
addMol("CuCl2|氯化銅", "Cu", "-", "-", "-", "620", "993", [{elem:"Cu",x:0,y:0,z:0,r:18,lpCount:0}, {elem:"Cl",x:-80,y:0,z:0,r:35,lpCount:0}, {elem:"Cl",x:80,y:0,z:0,r:35,lpCount:0}], [[0, 1, "ionic_thin"], [0, 2, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">🧂 物質性質</div><div class="info-body"><strong>氯化銅 (CuCl₂)</strong><br>燃燒呈藍綠色火焰。</div></div>');
addMol("PbI2|碘化鉛", "Pb", "-", "-", "-", "402", "953", [{elem:"Pb",x:0,y:0,z:0,r:25,lpCount:0}, {elem:"I",x:-90,y:0,z:0,r:40,lpCount:0}, {elem:"I",x:90,y:0,z:0,r:40,lpCount:0}], [[0, 1, "ionic_thin"], [0, 2, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">✨ 物質性質</div><div class="info-body"><strong>碘化鉛 (PbI₂)</strong><br>亮黃色晶體，用於「黃金雨」實驗。</div></div>');
addMol("CaH2|氫化鈣", "Ca", "-", "-", "-", "816", "分解", [{elem:"Ca",x:0,y:0,z:0,r:22,lpCount:0}, {elem:"H",x:-70,y:0,z:0,r:15,lpCount:0}, {elem:"H",x:70,y:0,z:0,r:15,lpCount:0}], [[0, 1, "ionic_thin"], [0, 2, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">⛺ 物質性質</div><div class="info-body"><strong>氫化鈣 (CaH₂)</strong><br>攜帶方便的氫氣發生劑。</div></div>');
addMol("Na2O|氧化鈉", "Na", "-", "-", "-", "1132", "1950", [{elem:"O",x:0,y:0,z:0,r:25,lpCount:0}, {elem:"Na",x:-90,y:0,z:0,r:20,lpCount:0}, {elem:"Na",x:90,y:0,z:0,r:20,lpCount:0}], [[0, 1, "ionic_thin"], [0, 2, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">🧂 物質性質</div><div class="info-body"><strong>氧化鈉 (Na₂O)</strong><br>鹼性氧化物。</div></div>');
addMol("K2O|氧化鉀", "K", "-", "-", "-", "740", "分解", [{elem:"O",x:0,y:0,z:0,r:25,lpCount:0}, {elem:"K",x:-100,y:0,z:0,r:22,lpCount:0}, {elem:"K",x:100,y:0,z:0,r:22,lpCount:0}], [[0, 1, "ionic_thin"], [0, 2, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">🧂 物質性質</div><div class="info-body"><strong>氧化鉀 (K₂O)</strong><br>極易吸濕，用於肥料計算基準。</div></div>');
addMol("Na2S|硫化鈉", "Na", "-", "-", "-", "1176", "-", [{elem:"S",x:0,y:0,z:0,r:30,lpCount:0}, {elem:"Na",x:-90,y:0,z:0,r:20,lpCount:0}, {elem:"Na",x:90,y:0,z:0,r:20,lpCount:0}], [[0, 1, "ionic_thin"], [0, 2, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">🧪 物質性質</div><div class="info-body"><strong>硫化鈉 (Na₂S)</strong><br>俗稱臭鹼，用於造紙與皮革工業。</div></div>');
addMol("FeCl2|氯化亞鐵", "Fe", "-", "-", "-", "677", "1023", [{elem:"Fe",x:0,y:0,z:0,r:18,lpCount:0}, {elem:"Cl",x:-85,y:0,z:0,r:35,lpCount:0}, {elem:"Cl",x:85,y:0,z:0,r:35,lpCount:0}], [[0, 1, "ionic_thin"], [0, 2, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">🧪 物質性質</div><div class="info-body"><strong>氯化亞鐵 (FeCl₂)</strong><br>淺綠色晶體，具還原性。</div></div>');
addMol("FeCl3|氯化鐵", "Fe", "-", "-", "-", "306", "315", [{elem:"Fe",x:0,y:0,z:0,r:18,lpCount:0}, {elem:"Cl",x:0,y:90,z:0,r:35,lpCount:0}, {elem:"Cl",x:-78,y:-45,z:0,r:35,lpCount:0}, {elem:"Cl",x:78,y:-45,z:0,r:35,lpCount:0}], [[0, 1, "ionic_thin"], [0, 2, "ionic_thin"], [0, 3, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">🧪 物質性質</div><div class="info-body"><strong>氯化鐵 (FeCl₃)</strong><br>黑棕色晶體，用於電路板蝕刻。</div></div>');
addMol("AlCl3|氯化鋁", "Al", "-", "-", "-", "192", "180 (昇華)", [{elem:"Al",x:0,y:0,z:0,r:18,lpCount:0}, {elem:"Cl",x:0,y:90,z:0,r:35,lpCount:0}, {elem:"Cl",x:-78,y:-45,z:0,r:35,lpCount:0}, {elem:"Cl",x:78,y:-45,z:0,r:35,lpCount:0}], [[0, 1, "ionic_thin"], [0, 2, "ionic_thin"], [0, 3, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">🧪 物質性質</div><div class="info-body"><strong>氯化鋁 (AlCl₃)</strong><br>路易斯酸，有機合成催化劑。</div></div>');
addMol("Al2O3|氧化鋁|剛玉", "Al", "-", "-", "-", "2072", "2977", [{elem:"O",x:-90,y:-20,z:0,r:25,lpCount:0}, {elem:"Al",x:-45,y:40,z:0,r:18,lpCount:0}, {elem:"O",x:0,y:-20,z:0,r:25,lpCount:0}, {elem:"Al",x:45,y:40,z:0,r:18,lpCount:0}, {elem:"O",x:90,y:-20,z:0,r:25,lpCount:0}], [[0, 1, "ionic_thin"], [1, 2, "ionic_thin"], [2, 3, "ionic_thin"], [3, 4, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">💎 物質性質</div><div class="info-body"><strong>氧化鋁 (Al₂O₃)</strong><br>剛玉，紅寶石與藍寶石的主要成分，硬度高。</div></div>');
addMol("Fe2O3|氧化鐵|赤鐵礦", "Fe", "-", "-", "-", "1565", "-", [{elem:"O",x:-90,y:-20,z:0,r:25,lpCount:0}, {elem:"Fe",x:-45,y:40,z:0,r:18,lpCount:0}, {elem:"O",x:0,y:-20,z:0,r:25,lpCount:0}, {elem:"Fe",x:45,y:40,z:0,r:18,lpCount:0}, {elem:"O",x:90,y:-20,z:0,r:25,lpCount:0}], [[0, 1, "ionic_thin"], [1, 2, "ionic_thin"], [2, 3, "ionic_thin"], [3, 4, "ionic_thin"]], null, '<div class="info-section"><div class="info-title">🧱 物質性質</div><div class="info-body"><strong>氧化鐵 (Fe₂O₃)</strong><br>紅棕色粉末，俗稱鐵鏽或紅土，為赤鐵礦成分。</div></div>');


// ==========================================
// 碳簇生成器 (放在檔案最後面)
// ==========================================
function addFullerene(name, n, radius, stretch = 1) {
    const atoms = [];
    if (n === 60) {
        const phi = (1 + Math.sqrt(5)) / 2; 
        const scale = 28;
        let rawVerts = [];
        const groups = [[0, 1, 3*phi], [1, 2+phi, 2*phi], [phi, 2, 2*phi+1]];
        groups.forEach(g => {
            const cycles = [[g[0],g[1],g[2]], [g[1],g[2],g[0]], [g[2],g[0],g[1]]];
            cycles.forEach(p => {
                for(let i=0; i<8; i++) {
                    let x = p[0] * ((i & 1) ? 1 : -1);
                    let y = p[1] * ((i & 2) ? 1 : -1);
                    let z = p[2] * ((i & 4) ? 1 : -1);
                    rawVerts.push({x, y, z});
                }
            });
        });
        const threshold = 0.1;
        rawVerts.forEach(v => {
            if(!atoms.some(a => Math.abs(a.x - v.x*scale) < threshold && Math.abs(a.y - v.y*scale) < threshold && Math.abs(a.z - v.z*scale) < threshold)) {
                atoms.push({elem: "C", x: v.x*scale, y: v.y*scale, z: v.z*scale, r: 5, lpCount: 0});
            }
        });
    } else {
        const phi = Math.PI * (3 - Math.sqrt(5)); 
        for(let i = 0; i < n; i++) {
            const y = 1 - (i / (n - 1)) * 2; 
            const r = Math.sqrt(1 - y * y);
            const theta = phi * i;
            const x = Math.cos(theta) * r;
            const z = Math.sin(theta) * r;
            atoms.push({elem: "C", x: x * radius, y: y * radius * stretch, z: z * radius, r: 5, lpCount: 0});
        }
    }
    const bonds = []; const bondSet = new Set();
    for(let i = 0; i < atoms.length; i++) {
        let dists = [];
        for(let j = 0; j < atoms.length; j++) {
            if(i === j) continue;
            let d = (atoms[i].x - atoms[j].x)**2 + (atoms[i].y - atoms[j].y)**2 + (atoms[i].z - atoms[j].z)**2;
            dists.push({id: j, d: d});
        }
        dists.sort((a, b) => a.d - b.d);
        for(let k = 0; k < 3; k++){
            let neighbor = dists[k].id;
            let id1 = Math.min(i, neighbor); let id2 = Math.max(i, neighbor);
            let key = `${id1}-${id2}`;
            if(!bondSet.has(key)){
                let type = (k === 0) ? "double" : "single"; 
                bonds.push([id1, id2, type]); bondSet.add(key);
            }
        }
    }
    // [修正] 補上 mp 和 bp 參數 ("-", "-")，避免參數錯位
    addMol(name, "C", "sp²", ["球狀", n===60?"Truncated Icosahedron":"Fullerene"], "120°", "-", "-", atoms, bonds);
}

// [修正] 加上中文名稱，讓「猜你想看」顯示更完整
addFullerene("C36|碳36", 36, 130);
addFullerene("C40|碳40", 40, 140);
addFullerene("C50|碳50", 50, 160);
addFullerene("C70|碳70", 70, 160, 1.4);
addFullerene("C80|碳80", 80, 170);
addFullerene("C100|碳100", 100, 190);

// Defensive cleanup for any "複合"/composite mode artifacts possibly defined here
(function(){
    const keys = ["模式: 複合","模式_複合","複合模式","compositeMode","modeComposite","composite"];
    if (typeof window !== 'undefined') {
        keys.forEach(k=>{
            try { if (window.hasOwnProperty(k)) delete window[k]; } catch(e){ try { window[k]=undefined; } catch(_){} }
        });
    }
})();




