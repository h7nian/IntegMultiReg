// Parse the actual TeX in generated HTML with KaTeX, without needing a browser.
const fs = require('fs');
const path = require('path');
const parse5 = require('parse5');
const katex = require('katex');
const root = path.resolve(process.argv[2] || 'docs');
let expressions = 0;
let failures = 0;
const skipped = new Set(['script', 'style', 'pre', 'code', 'textarea', 'noscript']);

function text(node) {
  return node.nodeName === '#text' ? node.value : (node.childNodes || []).map(text).join('');
}
function validate(tex, displayMode, file) {
  try {
    katex.renderToString(tex, {displayMode, throwOnError: true, trust: false});
    expressions++;
  } catch (error) {
    console.error(path.relative(root, file) + ': ' + error.message);
    failures++;
  }
}
function visit(node, file) {
  const attrs = (node.attrs || []).reduce((out, a) => {out[a.name] = a.value; return out;}, {});
  const classes = (attrs.class || '').split(/\s+/);
  if (classes.includes('katex')) return;
  if (classes.includes('math') || classes.includes('reqn')) {
    let tex = text(node).trim();
    tex = tex.replace(/^\\\(|\\\)$/g, '').replace(/^\\\[|\\\]$/g, '');
    validate(tex, classes.includes('display'), file);
    return;
  }
  if (skipped.has(node.tagName)) return;
  if (node.nodeName === '#text') {
    // Rd uses $$...$$ and \(...\); Pandoc can also emit \[...\].
    const pattern = /\$\$([\s\S]*?)\$\$|\\\(([\s\S]*?)\\\)|\\\[([\s\S]*?)\\\]/g;
    let match;
    while ((match = pattern.exec(node.value)) !== null) {
      validate(match[1] !== undefined ? match[1] : match[2] !== undefined ? match[2] : match[3],
               match[2] === undefined, file);
    }
  }
  (node.childNodes || []).forEach(child => visit(child, file));
}
function walk(dir) {
  fs.readdirSync(dir).forEach(name => {
    const file = path.join(dir, name);
    if (fs.statSync(file).isDirectory()) walk(file);
    else if (name.endsWith('.html')) {
      const before = expressions;
      visit(parse5.parse(fs.readFileSync(file, 'utf8')), file);
      const rd = path.join(__dirname, '..', '..', 'man', name.replace(/\.html$/, '.Rd'));
      if (fs.existsSync(rd)) {
        const expected = (fs.readFileSync(rd, 'utf8').match(/\\(?:eqn|deqn)\{/g) || []).length;
        if (expressions - before < expected) {
          console.error(path.relative(root, file) + ': missing rendered TeX (' +
                        (expressions - before) + ' found; ' + expected + ' documented).');
          failures++;
        }
      }
    }
  });
}
walk(root);
if (!expressions && !failures) {
  console.error('No mathematical expressions found in the generated website.');
  failures++;
}
if (failures) process.exit(1);
console.log('KaTeX parsed ' + expressions + ' mathematical expressions without errors.');
