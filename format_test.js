const fs = require('fs');

const files = [
  'c:/Jeffrey/Chem Proj/index.html',
  'c:/Jeffrey/Chem Proj/data_molecules.js'
];

let text = files.map(p => fs.readFileSync(p,'utf8')).join('\n');

// Find quoted tokens that end with + or - (including variants like ^2-)
const tokenRe = /["']([^"']+[+-])["']/g;
let m; const tokens = new Set();
while ((m = tokenRe.exec(text)) !== null) tokens.add(m[1]);

function parseAtoms(f) {
  let formula = f.replace(/\[/g,'(').replace(/\]/g,')').replace(/\{/g,'(').replace(/\}/g,')');
  const stack = [{}];
  let i = 0;
  while (i < formula.length) {
    let char = formula[i];
    if (char === '(') { stack.push({}); i++; }
    else if (char === ')') {
      i++; let mStr = "";
      while (i < formula.length && /\d/.test(formula[i])) { mStr += formula[i]; i++; }
      let m = parseInt(mStr || "1"), popped = stack.pop(), top = stack[stack.length - 1];
      for (let el in popped) { top[el] = (top[el] || 0) + popped[el] * m; }
    } else {
      let match = formula.slice(i).match(/^([A-Z][a-z]*)(\d*)/);
      if (match) {
        let el = match[1], count = parseInt(match[2] || "1"), top = stack[stack.length - 1];
        top[el] = (top[el] || 0) + count; i += match[0].length;
      } else { i++; }
    }
  }
  return stack[0];
}

function formatFormulaConsole(str) {
  if (!str) return "";
  let s = str.trim();
  let mainPart = s;
  let chargePart = "";
  s = s.replace(/\^/g, '').replace(/\s+/g, '');
  const m = s.match(/^(.*?)(\d+([+-])|[+-])$/);
  if (m) {
    let rawCharge = m[2] || m[0].slice(m[1].length);
    const startIdx = s.length - rawCharge.length;
    mainPart = s.slice(0, startIdx);
    const digitMatch = rawCharge.match(/^(\d+)/);
    if (digitMatch) {
      const digits = digitMatch[1];
      const charBeforeDigits = s[startIdx - 1];
      // If multiple leading digits and preceding char is a letter, assume
      // all but the last digit belong to the formula (e.g. O22- -> O2 + 2-)
      if (digits.length > 1 && charBeforeDigits && /[A-Za-z]/.test(charBeforeDigits)) {
        const leading = digits.slice(0, -1);
        const lastDigit = digits.slice(-1);
        mainPart += leading;
        rawCharge = lastDigit + rawCharge.slice(digits.length);
      } else {
        try {
          const elCounts = parseAtoms(mainPart);
          const distinctEls = Object.keys(elCounts).length;
          if (charBeforeDigits && /[A-Za-z]/.test(charBeforeDigits) && distinctEls > 1) {
            const signOnly = rawCharge.slice(digits.length);
            mainPart += digits;
            rawCharge = signOnly || '';
          }
        } catch (e) {}
      }
    }
    if (rawCharge) {
      const chargeDisplay = rawCharge.replace(/(\d+)/g, '$1').replace(/([+-]+)/g, '<charge>$1</charge>');
      chargePart = '^' + chargeDisplay + '^';
    }
  }
  mainPart = mainPart.replace(/(\d+)(?![^<]*>)/g, '_$1_');
  return mainPart + chargePart;
}

console.log('Found', tokens.size, 'candidate tokens');
for (let t of tokens) {
  console.log(t, '->', formatFormulaConsole(t));
}
