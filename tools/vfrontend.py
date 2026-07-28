"""Verify frontend/index.html landed intact. Run from D:\\dnainsight."""
import hashlib
import os
import re
import subprocess
import sys

TARGET = os.path.join('frontend', 'index.html')
APPJS = os.path.join('tools', 'app.js')
EXP_LINES = 2421
EXP_BYTES = 108961
EXP_MD5 = 'c9f59788377b6f353debee3f1dfd9ffc'

fails = []
infos = []


def check(ok, msg):
    print(('  [ok]   ' if ok else '  [FAIL] ') + msg)
    if not ok:
        fails.append(msg)
    return ok


raw = open(TARGET, 'rb').read()
text = raw.decode('utf-8')

print('=== 1. size and line count ===')
nlines = text.count('\n')
nbytes = len(raw)
check(nlines == EXP_LINES, 'line count is %d, expected %d' % (nlines, EXP_LINES))

if nbytes == EXP_BYTES:
    check(True, 'byte length is %d, exactly as expected' % nbytes)
else:
    delta = nbytes - EXP_BYTES
    crlf = (delta == EXP_LINES)
    note = 'consistent with CRLF conversion (every LF became CRLF)' if crlf \
        else 'NOT consistent with CRLF conversion'
    check(False, 'byte length is %d, expected %d, delta %+d, %s'
          % (nbytes, EXP_BYTES, delta, note))

print('  [info] CRLF present in file: %s' % ('yes' if b'\r\n' in raw else 'no'))
md5 = hashlib.md5(raw).hexdigest()
check(md5 == EXP_MD5, 'md5 is %s, expected %s' % (md5, EXP_MD5))

print('=== 2. document structure ===')
m = re.search(r'<script>(.*?)</script>', text, re.S)
if not m:
    check(False, 'could not find the inline <script> block')
    print('\nRESULT: FAIL')
    sys.exit(1)
js = m.group(1)
markup = text[:m.start()] + text[m.end():]
print('  [info] inline script is %d chars; real markup is %d chars' % (len(js), len(markup)))

for tag, pat in [('<html', r'<html'), ('</html>', r'</html>'),
                 ('<body>', r'<body>'), ('</body>', r'</body>')]:
    in_markup = len(re.findall(pat, markup))
    in_full = len(re.findall(pat, text))
    check(in_markup == 1,
          'exactly one %s in real markup (found %d; %d in whole file, '
          'the extras are inside JS strings and are correct)' % (tag, in_markup, in_full))

print('=== 3. node --check on the extracted script ===')
with open(APPJS, 'w', encoding='utf-8', newline='\n') as fh:
    fh.write(js)
try:
    p = subprocess.run(['node', '--check', APPJS], capture_output=True, text=True, shell=True)
    out = ((p.stdout or '') + (p.stderr or '')).strip()
    check(p.returncode == 0, 'node --check exit %d %s'
          % (p.returncode, ('output: ' + out) if out else '(no output, clean parse)'))
except Exception as exc:
    check(False, 'could not run node --check: %s' % exc)

print('=== 4. no browser storage in executable code ===')
stripped = re.sub(r'/\*.*?\*/', '', js, flags=re.S)
stripped = re.sub(r'(?<!:)//[^\n]*', '', stripped)
for word in ('localStorage', 'sessionStorage', 'indexedDB'):
    n_code = stripped.count(word)
    n_all = js.count(word)
    check(n_code == 0, '%s absent from executable code (%d in code, %d in the '
          'whole script including the comment that forbids it)' % (word, n_code, n_all))

print('=== 5. no em dash or en dash ===')
for name, ch in (('em dash', '\u2014'), ('en dash', '\u2013')):
    check(text.count(ch) == 0, 'zero %s characters (found %d)' % (name, text.count(ch)))

print('=== 6. inline handlers are all defined ===')
handlers = re.findall(r'\b(?:onclick|oninput|onchange)\s*=\s*"\s*([A-Za-z_$][A-Za-z0-9_$]*)\s*\(',
                      text)
uniq = sorted(set(handlers))
print('  [info] %d handler attributes, %d distinct names' % (len(handlers), len(uniq)))
missing = [n for n in uniq if ('function ' + n) not in text]
check(not missing, 'every handler name is defined as "function NAME"%s'
      % ('' if not missing else '; UNDEFINED: ' + ', '.join(missing)))

other = re.findall(r'\b(?:onclick|oninput|onchange)\s*=\s*"\s*(?![A-Za-z_$][A-Za-z0-9_$]*\s*\()([^"]{0,60})',
                   text)
for o in other:
    infos.append('handler attribute not a plain function call, skipped: %s' % o.strip())

print('=== 7. required contract strings ===')
required = ['goodRepute', 'badRepute', 'noRepute', '#60B060', '#FF9090', '#C0C0C0',
            '#998EC3', '#F1A340', 'SLIDER_DEFS', 'CLINVAR_CODES', 'min_stars',
            'facet-genes', 'facet-topics', 'facet-medicines', 'facet-conditions',
            'cnt-visible', 'doubleAllowed', 'pop-select', 'onPopulationChange',
            'renderQcBanner', 'renderGenosets', 'renderPgx', 'renderTraits',
            'renderPrs', 'renderQc', 'SNPs only', 'resetFilters', 'magnitude_factors',
            'Both calls are kept', 'Not testable on your array', 'CC BY-NC-SA 3.0 US',
            'findings/v2', 'aria-label']
absent = [s for s in required if s not in text]
check(not absent, 'all %d required strings present%s'
      % (len(required), '' if not absent else '; MISSING: ' + ', '.join(absent)))

print()
for i in infos:
    print('  [info] ' + i)
print()
if fails:
    print('RESULT: FAIL, %d problem(s):' % len(fails))
    for f in fails:
        print('  - ' + f)
    sys.exit(1)
print('RESULT: PASS, all checks green')
