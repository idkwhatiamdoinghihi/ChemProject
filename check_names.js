const fs = require('fs');
const path = 'c:/Jeffrey/Chem Proj/index.html';
const txt = fs.readFileSync(path,'utf8');

// Extract COMMON_EN_NAMES block
const m = txt.match(/const\s+COMMON_EN_NAMES\s*=\s*\{([\s\S]*?)\};/);
let common = {};
if (m) {
  const body = m[1];
  const kvRe = /['"]([^'"]+)['"]\s*:\s*['"]([^'"]+)['"]/g;
  let mm;
  while ((mm = kvRe.exec(body)) !== null) {
    common[mm[1].toUpperCase()] = mm[2];
  }
}

// Find candidate keys (strings containing +, -, or ^)
const keyRe = /["']([^"']*[+\-\^][^"']*)["']\s*:/g;
let kk; const keys = new Set();
while ((kk = keyRe.exec(txt)) !== null) {
  keys.add(kk[1]);
}

function norm(k){ return String(k||'').replace(/\^/g,'').replace(/\s+/g,'').toUpperCase(); }

function deriveName(key){
  const parts = String(key).split('|').map(s=>s.trim()).filter(Boolean);
  let nameEN = '';
  for (let i=1;i<parts.length;i++){
    if (/[A-Za-z]/.test(parts[i])) { nameEN = parts[i]; break; }
  }
  if (!nameEN){
    const p0 = parts[0]||''; const p1 = parts[1]||parts[0]||'';
    const tryKeys = [p1, p0, norm(p0), norm(p1), p0.replace(/\s+/g,''), p0.replace(/\^/g,'')];
    for (let tk of tryKeys){ if (!tk) continue; if (common[tk.toUpperCase()]) { nameEN = common[tk.toUpperCase()]; break; } }
  }
  if (!nameEN) nameEN = parts[0]||key;
  return nameEN;
}

const results = [];
keys.forEach(k=>{
  results.push({key:k, name: deriveName(k), mapped: !!common[norm(k)]});
});

results.sort((a,b)=>a.key.localeCompare(b.key));
console.log('Found', results.length, 'ion/formula-like keys');
results.forEach(r=> console.log(r.key, '=>', r.name, r.mapped ? '(mapped)' : ''));
