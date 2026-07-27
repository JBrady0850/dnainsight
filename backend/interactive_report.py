"""
interactive_report.py -- the self-contained offline interactive report.

WHAT THIS IS
------------
One HTML file containing both the findings AND a working filter engine. It opens
with no server running, no internet, and no DNAInsight installation. Give it to
yourself on a USB stick in ten years and it still works.

This is the equivalent of the reference product's single-file interactive output,
and it is the artifact a user actually keeps. The two static reports
(genetic_report.py, doctor_report.py) are for printing and for handing to a
clinician; this one is for exploring.

DESIGN CONSTRAINTS, ALL DELIBERATE
----------------------------------
1. NO external requests. No CDN, no web fonts, no analytics. The whole point is
   that it works offline and leaks nothing. Everything is inlined, and the chart
   is hand-drawn SVG rather than a charting library.
2. NO browser storage. Same project-wide rule as the main app.
3. The data is embedded as a JSON island in a <script type="application/json">
   tag, not interpolated into JavaScript source. That keeps a stray apostrophe in
   an interpretation string from breaking the whole file, and it means the data
   can be lifted back out with a parser.
4. Filtering here is CLIENT side, which is the opposite of the main app. That is
   correct: there is no server to ask, and the export is a fixed snapshot.
5. Every caveat, every strand warning and every "not testable" distinction
   carries through. An offline copy must not be more confident than the live app.

SECURITY NOTE
-------------
Every value that reaches the document goes through _esc(). A raw DNA file is
attacker-controllable in principle (someone hands you a file), so an
interpretation string is never trusted as markup. The JSON island additionally
escapes "<" so a crafted string cannot close the script tag early.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from typing import Any, Iterable

from . import APP_VERSION

__all__ = ["generate_interactive_report", "build_payload"]


# ---------------------------------------------------------------------------
# Escaping
# ---------------------------------------------------------------------------

_HTML_ESCAPES = {
    "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;",
}


def _esc(value: Any) -> str:
    """HTML-escape any value. Never trust a string that came from a data file."""
    if value is None:
        return ""
    return "".join(_HTML_ESCAPES.get(c, c) for c in str(value))


def _json_island(data: Any) -> str:
    """Serialise for embedding in a script tag, safely.

    Escaping "<" prevents a crafted value containing "</script>" from closing
    the tag early, which would break the document and could inject markup.
    """
    text = json.dumps(data, ensure_ascii=False, separators=(",", ":"))
    return text.replace("<", "\\u003c").replace("\u2028", "\\u2028") \
               .replace("\u2029", "\\u2029")


# ---------------------------------------------------------------------------
# Payload
# ---------------------------------------------------------------------------

# Only these keys travel into the file. Trimming matters: a full scan can carry
# tens of thousands of findings, and population_series alone is 16 objects each.
_KEEP = (
    "rsid", "entity_type", "gene", "chromosome", "position", "genotype",
    "token", "zygosity", "magnitude", "magnitude_source", "magnitude_factors",
    "repute", "summary", "interpretation", "confidence", "clinical_sig",
    "clinvar_sig_code", "review_stars", "cpic_level", "evidence",
    "publications", "conditions", "conditions_list", "silo", "category",
    "freq", "freq_population", "freq_band", "freq_color", "freq_derived",
    "freq_method", "freq_flipped", "freq_ambiguous", "gmaf", "minor_allele",
    "population_series", "orientation", "stabilized_orientation", "flipped",
    "ambiguous", "dubious", "variant_allele", "variant_copies", "carrier",
    "count", "labels", "conflict", "calls", "comparison", "topics",
    "medicines", "criteria", "matched_rsids", "coverage", "percentile",
    "band", "reliable", "caveats", "partial_coverage", "name",
)


def build_payload(profile: dict, findings: Iterable[dict],
                  extras: dict | None = None) -> dict:
    """Assemble the JSON payload that gets embedded in the report."""
    extras = extras or {}
    trimmed = []
    for f in findings or []:
        row = {k: f.get(k) for k in _KEEP if k in f}
        trimmed.append(row)

    return {
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "app_version": APP_VERSION,
        "patient": {
            "name": profile.get("name", ""),
            "dob": profile.get("dob", ""),
            "sex": profile.get("sex", ""),
            "provider": profile.get("provider", ""),
        },
        "population": extras.get("population", ""),
        "counts": extras.get("counts", {}),
        "qc": extras.get("qc", {}),
        "sources": extras.get("sources", []),
        "genoset_incomplete": [
            {"rsid": g.get("rsid"), "aka": g.get("aka", ""),
             "summary": g.get("summary", ""), "coverage": g.get("coverage")}
            for g in (extras.get("genosets", {}) or {}).get("incomplete", [])
        ],
        "blood_type": extras.get("blood_type", {}),
        "findings": trimmed,
    }


# ---------------------------------------------------------------------------
# The document
# ---------------------------------------------------------------------------

_CSS = """
:root{--good:#60B060;--bad:#FF9090;--unset:#C0C0C0;--blue:#1a3a6b;
--mid:#2980b9;--bg:#f4f6f9;--text:#2c3e50;--muted:#7f8c8d;--line:#e3e9ef}
body.cb{--good:#998EC3;--bad:#F1A340}
*{box-sizing:border-box;margin:0;padding:0}
body{font-family:'Segoe UI',Arial,sans-serif;background:var(--bg);
color:var(--text);font-size:15px;line-height:1.45}
header{background:var(--blue);color:#fff;padding:14px 20px}
header h1{font-size:1.2em}
header .meta{font-size:.8em;opacity:.85;margin-top:3px}
.wrap{display:flex;align-items:flex-start;gap:16px;padding:16px 20px}
.panel{width:270px;flex-shrink:0;background:#fff;border-radius:8px;padding:14px;
box-shadow:0 2px 8px rgba(0,0,0,.1);position:sticky;top:14px;max-height:92vh;
overflow-y:auto}
.panel h3{font-size:.7em;text-transform:uppercase;letter-spacing:1px;
color:var(--muted);margin:14px 0 6px;border-top:1px solid var(--line);padding-top:10px}
.panel h3:first-of-type{border:none;margin-top:0;padding-top:0}
.main{flex:1;min-width:0}
.card{background:#fff;border-radius:8px;padding:14px 16px;margin-bottom:12px;
box-shadow:0 2px 8px rgba(0,0,0,.1)}
.geno{background:#fff;border-radius:8px;border-left:4px solid var(--unset);
padding:11px 14px;margin-bottom:9px;box-shadow:0 1px 5px rgba(0,0,0,.09);
display:grid;grid-template-columns:minmax(0,1fr) 230px;gap:14px}
.geno.g{border-left-color:var(--good);background:#f2fbf2}
.geno.b{border-left-color:var(--bad);background:#fff5f5}
.geno.dub{border-left-style:dashed}
@media(max-width:900px){.geno{grid-template-columns:1fr}.wrap{flex-direction:column}
.panel{width:100%;position:static;max-height:none}}
.gh{display:flex;gap:8px;align-items:center;flex-wrap:wrap;margin-bottom:5px}
.gid{font-family:Consolas,monospace;font-weight:700}
.tok{font-family:Consolas,monospace;color:var(--mid)}
.pill{background:var(--blue);color:#fff;border-radius:11px;padding:1px 8px;
font-size:.78em;font-weight:700}
.pill.z{background:#b9c1c9}.pill.l{background:#8fa6bd}
.pill.m{background:var(--mid)}.pill.h{background:#b5341f}
.bdg{font-size:.71em;padding:2px 7px;border-radius:9px;background:#ecf0f1;color:#556}
.bdg.rx{background:#fdecea;color:#8b1c13}.bdg.ac{background:#fff4e5;color:#8a4b08}
.bdg.wn{background:#fff8e1;color:#7d5900}.bdg.cp{background:#e8f0fe;color:#17408b}
.bdg.st{background:#fef6d8;color:#7a5c00}
.sum{font-weight:600;margin-bottom:3px}
.body{font-size:.87em;color:#445}
.crit{font-family:Consolas,monospace;font-size:.79em;background:#f7f9fb;
border:1px solid var(--line);border-radius:4px;padding:6px 8px;margin-top:7px;
white-space:pre-wrap;word-break:break-word}
table.kv{width:100%;font-size:.78em}
table.kv td{padding:2px 3px}
table.kv td:first-child{color:var(--muted);width:47%}
table.kv td:last-child{font-weight:600;text-align:right}
.bars{display:flex;gap:2px;align-items:flex-end;height:30px;margin-top:7px}
.bar{flex:1;min-height:2px;border-radius:2px 2px 0 0;background:#EFEFEF}
.bar.you{outline:2px solid var(--blue)}
.blab{font-size:.58em;color:var(--muted);text-align:center;margin-top:2px}
.tags{margin-top:7px;font-size:.77em;display:flex;flex-wrap:wrap;gap:4px;
align-items:center}
.tags .tl{color:var(--muted)}
.tag{background:#e8f4fd;color:var(--blue);border-radius:9px;padding:1px 7px}
.cav{margin-top:8px;background:#f7f9fb;border-left:3px solid var(--muted);
padding:6px 9px;font-size:.81em}
.cav ul{margin:3px 0 0 15px}
.cfl{margin-top:8px;background:#fff8e1;border:1px solid #ffc107;border-radius:4px;
padding:7px 9px;font-size:.81em}
.ban{padding:9px 13px;border-radius:6px;font-size:.86em;margin-bottom:12px}
.ban.i{background:#e8f4fd;color:var(--blue);border:1px solid #b8daf5}
.ban.w{background:#fff8e1;color:#7d5900;border:1px solid #ffc107}
label{font-size:.85em;display:block;margin-bottom:3px;cursor:pointer}
label.inl{display:flex;gap:5px;align-items:center}
input[type=text],select{width:100%;padding:5px 8px;border:1px solid #ccd;
border-radius:4px;font-family:inherit;font-size:.9em}
input[type=range]{width:100%}
.btn{padding:5px 11px;border-radius:5px;border:1px solid var(--mid);
background:#fff;color:var(--mid);cursor:pointer;font-family:inherit;font-size:.84em}
.btn:hover{background:#e8f4fd}
.cnt{display:flex;gap:6px;font-size:.78em;text-align:center;margin-top:6px}
.cnt div{flex:1;background:var(--bg);border-radius:4px;padding:4px}
.cnt strong{display:block;font-size:1.2em}
.muted{color:var(--muted)}
.mono{font-family:Consolas,monospace}
details summary{cursor:pointer}
.empty{text-align:center;padding:34px;color:var(--muted)}
footer{padding:16px 20px;font-size:.8em;color:var(--muted)}
@media print{.panel{display:none}.wrap{display:block}.geno{break-inside:avoid;
box-shadow:none;border:1px solid var(--line)}}
"""

_JS = r"""
/* Offline report engine. Client side by necessity: there is no server here.
   No external requests, no browser storage, no dependencies. */
var DATA = JSON.parse(document.getElementById('payload').textContent);
var F = DATA.findings || [];
var RAMP = [[50,'#FFFFFF'],[20,'#FFF3EE'],[10,'#FFE0D4'],[5,'#FFC5B2'],
            [1,'#FA9A80'],[0.1,'#EF6A4C'],[0,'#D8351B']];

function esc(v){if(v===null||v===undefined)return '';
  return String(v).replace(/[&<>"']/g,function(c){
    return {'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c];});}
function color(f){if(f===null||f===undefined||isNaN(f))return '#EFEFEF';
  for(var i=0;i<RAMP.length;i++){if(f>=RAMP[i][0])return RAMP[i][1];}
  return '#D8351B';}
/* A null magnitude counts as 1: unscored means nobody assessed it, which
   belongs above "assessed and boring" and below "interesting". */
function mag(f){var m=f.magnitude;return (m===null||m===undefined||isNaN(m))?1:Number(m);}
function magCls(m){if(m===null)return 'l';if(m===0)return 'z';
  if(m<2)return 'l';if(m<6)return 'm';return 'h';}
var EXEMPT={genoset:1,trait:1,prs:1};

var ST={q:'',minMag:0,maxMag:10,minPub:0,maxPub:0,minFreq:0,maxFreq:100,
        good:1,unset:1,bad:1,snp:1,genoset:1,trait:1,prs:1,
        gene:'',sort:'magnitude',order:'desc',allowed:50,requireFreq:0,
        ambigOnly:0,conflictOnly:0};

function bounds(){
  var mags=[],pubs=[],frq=[];
  F.forEach(function(f){
    mags.push(mag(f));
    if(typeof f.publications==='number')pubs.push(f.publications);
    if(typeof f.freq==='number')frq.push(f.freq);});
  return {mag:[0,Math.max.apply(null,mags.concat([10]))],
          pub:[0,pubs.length?Math.max.apply(null,pubs):0],
          frq:[0,frq.length?Math.max.apply(null,frq):100]};
}
var B=bounds();

function match(f){
  var ex=!!EXEMPT[f.entity_type];
  if(!ST[f.entity_type||'snp'])return false;
  var m=mag(f);
  if(m<ST.minMag||m>ST.maxMag)return false;
  var rep=f.repute==='Good'?'good':(f.repute==='Bad'?'bad':'unset');
  if(!ST[rep])return false;
  /* Frequency and publication filters do NOT apply to genosets, traits or
     scores. They have no single position, so they have no frequency. */
  if(!ex){
    var p=typeof f.publications==='number'?f.publications:0;
    if(p<ST.minPub||(ST.maxPub&&p>ST.maxPub))return false;
    if(ST.requireFreq&&(f.freq===null||f.freq===undefined))return false;
    if(typeof f.freq==='number'){
      if(f.freq<ST.minFreq||f.freq>ST.maxFreq)return false;
    }else if(ST.minFreq>0)return false;
  }
  if(ST.gene&&String(f.gene||'').toUpperCase()!==ST.gene)return false;
  if(ST.ambigOnly&&!(f.ambiguous||f.freq_ambiguous||f.flipped))return false;
  if(ST.conflictOnly&&!f.conflict)return false;
  if(ST.q){
    var hay=[f.rsid,f.gene,f.summary,f.interpretation,f.conditions,f.token,
             f.genotype,f.evidence,f.criteria,f.name]
      .concat(f.topics||[],f.medicines||[],f.conditions_list||[]).join(' ');
    var rx;
    try{rx=new RegExp(ST.q,'i');}
    catch(e){rx=new RegExp(ST.q.replace(/[.*+?^${}()|[\]\\]/g,'\\$&'),'i');}
    if(!rx.test(hay))return false;
  }
  return true;
}

function sorted(rows){
  var k=ST.sort,dir=ST.order==='asc'?1:-1;
  return rows.slice().sort(function(a,b){
    var x,y;
    if(k==='magnitude'){x=mag(a);y=mag(b);}
    else if(k==='frequency'){x=a.freq===null||a.freq===undefined?-1:a.freq;
      y=b.freq===null||b.freq===undefined?-1:b.freq;}
    else if(k==='publications'){x=a.publications||0;y=b.publications||0;}
    else if(k==='stars'){x=a.review_stars||0;y=b.review_stars||0;}
    else if(k==='gene'){x=String(a.gene||'~');y=String(b.gene||'~');}
    else {x=String(a.rsid||'');y=String(b.rsid||'');}
    if(x<y)return -1*dir; if(x>y)return 1*dir; return 0;});
}

function card(f){
  var m=f.magnitude===null||f.magnitude===undefined?null:Number(f.magnitude);
  var cls=f.repute==='Good'?'g':(f.repute==='Bad'?'b':'');
  var bd=[];
  if(f.silo==='pre_prescription')bd.push(['rx','prescription critical']);
  else if(f.silo==='actionable')bd.push(['ac','actionable']);
  if(f.entity_type!=='snp')bd.push(['',''+f.entity_type]);
  if(f.zygosity&&f.zygosity!=='no_call')bd.push(['',f.zygosity]);
  if(f.zygosity==='no_call')bd.push(['wn','no call']);
  if(f.review_stars)bd.push(['st',f.review_stars+' star'+(f.review_stars>1?'s':'')]);
  if(f.cpic_level)bd.push(['cp','CPIC '+f.cpic_level]);
  if(f.carrier===false)bd.push(['','not a carrier']);
  if(f.conflict)bd.push(['wn','files disagree']);
  if(f.flipped)bd.push(['','strand flipped']);
  if(f.ambiguous||f.freq_ambiguous)bd.push(['wn','strand ambiguous']);
  if(f.partial_coverage)bd.push(['wn','partial coverage']);

  var h='<div class="gh"><span class="gid">'+esc(f.rsid)+'</span>'+
    (f.token?'<span class="tok">'+esc(f.token)+'</span>':'')+
    (f.gene?'<span class="muted">'+esc(f.gene)+'</span>':'')+
    '<span class="pill '+magCls(m)+'">'+(m===null?'unscored':m)+'</span>'+
    bd.map(function(x){return '<span class="bdg '+x[0]+'">'+esc(x[1])+'</span>';}).join('')+
    '</div>';

  var kv=[],add=function(k,v){if(v===null||v===undefined||v===''||v===false)return;
    kv.push('<tr><td>'+esc(k)+'</td><td>'+v+'</td></tr>');};
  add('Magnitude',m===null?'unscored':m);
  add('Repute',f.repute||((f.entity_type==='trait'||f.entity_type==='prs')?
    '<span class="muted">neutral by design</span>':'<span class="muted">not set</span>'));
  add('Confidence',f.confidence&&f.confidence!=='none'?f.confidence:
    '<span class="muted">none</span>');
  if(f.entity_type==='snp'){
    add('Frequency',(f.freq===null||f.freq===undefined)?
      '<span class="muted">no data</span>':
      (Math.round(f.freq*100)/100)+'%'+(f.freq_derived?' <span class="muted">(derived)</span>':''));
    add('Band',f.freq_band&&f.freq_band!=='unknown'?f.freq_band.replace('_',' '):null);
    add('GMAF',f.gmaf===null||f.gmaf===undefined?null:f.gmaf);
    add('Publications',f.publications||null);
    add('Chromosome',f.chromosome||null);
    add('Position',f.position||null);
  }
  add('Evidence',f.evidence||null);
  add('Copies',(typeof f.variant_copies==='number')?f.variant_copies+' of 2':null);
  add('Coverage',(f.coverage!==null&&f.coverage!==undefined)?
    Math.round(f.coverage*100)+'%':null);
  add('Percentile',(f.percentile!==null&&f.percentile!==undefined)?f.percentile:null);
  add('Seen',f.count>1?f.count+' files':null);
  add('Strand',f.ambiguous?'cannot be verified':(f.flipped?'complemented':null));

  var bars='';
  if((f.population_series||[]).length){
    var mx=1;f.population_series.forEach(function(s){
      if(typeof s.frequency==='number'&&s.frequency>mx)mx=s.frequency;});
    bars='<div class="bars">'+f.population_series.map(function(s){
      var v=s.frequency,hh=(v===null||v===undefined)?2:Math.max(2,v/mx*28);
      var t=(v===null||v===undefined)?('No data for '+s.code):
        (Math.round(v*10)/10)+'% of '+s.code+(s.yours?' share your genotype':'');
      return '<div class="bar'+(s.yours?' you':'')+'" style="height:'+hh+
        'px;background:'+color(v)+'" title="'+esc(t)+'"></div>';}).join('')+
      '</div><div class="blab">'+esc(f.population_series.map(function(s){
        return s.code;}).join(' '))+'</div>';
  }

  var tags='';
  [['Topics',f.topics],['Medicines',f.medicines],
   ['Conditions',f.conditions_list]].forEach(function(p){
    if(p[1]&&p[1].length)tags+='<div class="tags"><span class="tl">'+p[0]+
      '</span>'+p[1].slice(0,12).map(function(v){
        return '<span class="tag">'+esc(v)+'</span>';}).join('')+'</div>';});

  var cav=[];
  (f.caveats||[]).forEach(function(c){cav.push(c);});
  if(f.ambiguous||f.freq_ambiguous)cav.push('This is an A/T or C/G genotype, so '+
    'the strand cannot be verified from the data. The reading shown may be the '+
    'complement.');
  if(m!==null&&m>=6&&f.entity_type==='snp')cav.push('Rare high impact calls on a '+
    'consumer array are sometimes miscalls. Confirm with a clinically validated '+
    'test before acting.');
  var cavHtml=cav.length?'<div class="cav"><strong>Read with care</strong><ul>'+
    cav.map(function(c){return '<li>'+esc(c)+'</li>';}).join('')+'</ul></div>':'';

  var cfl='';
  if(f.conflict&&(f.calls||[]).length){
    cfl='<div class="cfl"><strong>Your files disagree here.</strong> '+
      'Both calls are kept and neither was chosen.<table style="margin-top:5px">'+
      f.calls.map(function(c){return '<tr><td>'+esc(c.label)+'</td><td class="mono">'+
        esc(c.genotype)+'</td></tr>';}).join('')+'</table></div>';
  }

  var fact=(f.magnitude_factors||[]).length?
    '<details style="margin-top:6px;font-size:.79em"><summary class="muted">'+
    'How this magnitude was calculated</summary><ul style="margin:4px 0 0 15px">'+
    f.magnitude_factors.map(function(x){return '<li>'+esc(x)+'</li>';}).join('')+
    '</ul></details>':'';

  return '<div class="geno '+cls+(f.dubious?' dub':'')+'"><div>'+h+
    (f.summary?'<div class="sum">'+esc(f.summary)+'</div>':'')+
    (f.interpretation&&f.interpretation!==f.summary?
      '<div class="body">'+esc(f.interpretation)+'</div>':'')+
    (f.criteria?'<div class="crit">'+esc(f.criteria)+'</div>':'')+
    cfl+cavHtml+tags+fact+
    '</div><div><table class="kv">'+kv.join('')+'</table>'+bars+'</div></div>';
}

function render(){
  var rows=sorted(F.filter(match));
  var shown=rows.slice(0,ST.allowed);
  document.getElementById('vis').textContent=shown.length;
  document.getElementById('off').textContent=Math.max(0,rows.length-shown.length);
  document.getElementById('tot').textContent=rows.length;
  var g={Good:0,Bad:0,unset:0};
  rows.forEach(function(f){
    g[f.repute==='Good'?'Good':(f.repute==='Bad'?'Bad':'unset')]++;});
  document.getElementById('rep').innerHTML=
    '<span style="color:var(--good)">&#9679;</span> '+g.Good+' good &nbsp; '+
    '<span style="color:var(--unset)">&#9679;</span> '+g.unset+' not set &nbsp; '+
    '<span style="color:var(--bad)">&#9679;</span> '+g.Bad+' bad';
  document.getElementById('list').innerHTML=shown.length?
    shown.map(card).join(''):
    '<div class="empty"><strong>Nothing matches these filters.</strong><br>'+
    '<button class="btn" style="margin-top:9px" onclick="reset()">Reset</button></div>';
  document.getElementById('more').innerHTML=(rows.length>shown.length)?
    '<button class="btn" onclick="ST.allowed*=2;render()">Show twice as many</button>':'';
}

function bind(){
  var g={};
  F.forEach(function(f){if(f.gene)g[f.gene]=(g[f.gene]||0)+1;});
  var opts=Object.keys(g).sort().map(function(k){
    return '<option value="'+esc(k.toUpperCase())+'">'+esc(k)+' ('+g[k]+')</option>';});
  document.getElementById('gene').innerHTML=
    '<option value="">All genes</option>'+opts.join('');
  ST.maxPub=B.pub[1];ST.maxFreq=B.frq[1];ST.maxMag=B.mag[1];
  ['mag','pub','frq'].forEach(function(k){
    var lo=document.getElementById(k+'lo'),hi=document.getElementById(k+'hi');
    var bb=B[k];
    lo.min=hi.min=bb[0];lo.max=hi.max=bb[1];lo.value=bb[0];hi.value=bb[1];
  });
  sync();
}

function sync(){
  ST.q=document.getElementById('q').value.trim();
  ST.minMag=+document.getElementById('maglo').value;
  ST.maxMag=+document.getElementById('maghi').value;
  if(ST.minMag>ST.maxMag){var t=ST.minMag;ST.minMag=ST.maxMag;ST.maxMag=t;}
  ST.minPub=+document.getElementById('publo').value;
  ST.maxPub=+document.getElementById('pubhi').value;
  ST.minFreq=+document.getElementById('frqlo').value;
  ST.maxFreq=+document.getElementById('frqhi').value;
  ['good','unset','bad','snp','genoset','trait','prs'].forEach(function(k){
    var e=document.getElementById('c_'+k); if(e)ST[k]=e.checked?1:0;});
  ST.requireFreq=document.getElementById('c_reqfreq').checked?1:0;
  ST.ambigOnly=document.getElementById('c_ambig').checked?1:0;
  ST.conflictOnly=document.getElementById('c_conf').checked?1:0;
  ST.gene=document.getElementById('gene').value;
  var s=document.getElementById('sort').value.split('|');
  ST.sort=s[0];ST.order=s[1];
  document.getElementById('magv').textContent=ST.minMag+' to '+ST.maxMag;
  document.getElementById('pubv').textContent=ST.minPub+' to '+ST.maxPub;
  document.getElementById('frqv').textContent=ST.minFreq+' to '+ST.maxFreq;
  render();
}

function reset(){
  document.getElementById('q').value='';
  ['c_good','c_unset','c_bad','c_snp','c_genoset','c_trait','c_prs'].forEach(
    function(id){var e=document.getElementById(id);if(e)e.checked=true;});
  ['c_reqfreq','c_ambig','c_conf'].forEach(function(id){
    var e=document.getElementById(id);if(e)e.checked=false;});
  document.getElementById('gene').value='';
  document.getElementById('sort').value='magnitude|desc';
  ST.allowed=50;
  bind();
}

function cb(on){document.body.classList.toggle('cb',on);render();}

document.addEventListener('keydown',function(e){
  if(e.key==='Escape')reset();});
window.addEventListener('DOMContentLoaded',bind);
"""


def generate_interactive_report(profile: dict, findings: Iterable[dict],
                                extras: dict | None = None) -> str:
    """Return the complete self-contained interactive report as one HTML string."""
    extras = extras or {}
    payload = build_payload(profile, findings, extras)
    qc = payload["qc"] or {}
    name = _esc(payload["patient"]["name"] or "Unnamed profile")

    qc_bits = []
    if qc.get("flipped"):
        qc_bits.append(f"{qc['flipped']} calls had their alleles complemented to "
                       "match the reference, which is routine")
    if qc.get("ambiguous"):
        qc_bits.append(f"{qc['ambiguous']} are A/T or C/G genotypes whose strand "
                       "cannot be verified from the data")
    if qc.get("no_call"):
        qc_bits.append(f"{qc['no_call']} were no calls and score zero")
    if qc.get("conflicts"):
        qc_bits.append(f"{qc['conflicts']} positions where two of your own files "
                       "disagree, with both readings kept")
    qc_html = ""
    if qc_bits:
        qc_html = ('<div class="ban w"><strong>Data quality.</strong> '
                   + _esc("; ".join(qc_bits)) + ".</div>")

    incomplete = payload["genoset_incomplete"]
    inc_html = ""
    if incomplete:
        rows = "".join(
            f'<tr><td class="mono">{_esc(g["rsid"])}</td>'
            f'<td>{_esc(g.get("aka") or "")}</td>'
            f'<td>{_esc(g.get("summary") or "")}</td></tr>'
            for g in incomplete[:200])
        inc_html = (
            '<div class="card"><h3 style="margin-bottom:6px">Rules that could not '
            f'be evaluated ({len(incomplete)})</h3>'
            '<div class="ban w">Each of these needs a position your array did not '
            'read, so it could not be evaluated at all. That is different from '
            'having been checked and found absent.</div>'
            '<table style="font-size:.84em"><tbody>' + rows + '</tbody></table></div>')

    slider = lambda key, label, note: (
        f'<h3>{label} <span class="muted" style="font-size:.85em">'
        f'<span id="{key}v"></span></span></h3>'
        + (f'<div class="muted" style="font-size:.74em;margin-bottom:3px">{note}</div>'
           if note else '')
        + f'<input type="range" id="{key}lo" oninput="sync()">'
        + f'<input type="range" id="{key}hi" oninput="sync()">')

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>DNAInsight interactive report, {name}</title>
<style>{_CSS}</style>
</head>
<body>
<header>
  <h1>DNAInsight interactive report</h1>
  <div class="meta">
    {name}
    {(" &middot; born " + _esc(payload["patient"]["dob"])) if payload["patient"]["dob"] else ""}
    {(" &middot; " + _esc(payload["patient"]["provider"])) if payload["patient"]["provider"] else ""}
    &middot; generated {_esc(payload["generated_at"])}
    &middot; DNAInsight v{_esc(payload["app_version"])}
    {(" &middot; frequencies for " + _esc(payload["population"])) if payload["population"] else ""}
  </div>
</header>

<div class="wrap">
  <aside class="panel">
    <h3>Search</h3>
    <input type="text" id="q" placeholder="gene, rsID, keyword" oninput="sync()">
    <div class="muted" style="font-size:.74em;margin-top:3px">
      Accepts a regular expression. Press Escape to reset everything.
    </div>

    {slider('mag', 'Magnitude', '')}
    {slider('pub', 'Publications', 'SNPs only, by design')}
    {slider('frq', 'Frequency %', 'SNPs only, by design')}

    <h3>Repute</h3>
    <label class="inl"><input type="checkbox" id="c_good" checked onchange="sync()"> Good</label>
    <label class="inl"><input type="checkbox" id="c_unset" checked onchange="sync()"> Not set</label>
    <label class="inl"><input type="checkbox" id="c_bad" checked onchange="sync()"> Bad</label>

    <h3>Show</h3>
    <label class="inl"><input type="checkbox" id="c_snp" checked onchange="sync()"> SNPs</label>
    <label class="inl"><input type="checkbox" id="c_genoset" checked onchange="sync()"> Genosets</label>
    <label class="inl"><input type="checkbox" id="c_trait" checked onchange="sync()"> Traits</label>
    <label class="inl"><input type="checkbox" id="c_prs" checked onchange="sync()"> Risk scores</label>

    <h3>Restrict</h3>
    <label class="inl"><input type="checkbox" id="c_reqfreq" onchange="sync()"> Has a frequency</label>
    <label class="inl"><input type="checkbox" id="c_ambig" onchange="sync()"> Strand flagged only</label>
    <label class="inl"><input type="checkbox" id="c_conf" onchange="sync()"> File disagreements only</label>

    <h3>Gene</h3>
    <select id="gene" onchange="sync()"></select>

    <h3>Sort</h3>
    <select id="sort" onchange="sync()">
      <option value="magnitude|desc">Magnitude, high to low</option>
      <option value="magnitude|asc">Magnitude, low to high</option>
      <option value="frequency|desc">Frequency, common first</option>
      <option value="frequency|asc">Frequency, rare first</option>
      <option value="publications|desc">Publications, most first</option>
      <option value="stars|desc">Review stars, most first</option>
      <option value="gene|asc">Gene, A to Z</option>
      <option value="rsid|asc">Identifier, ascending</option>
    </select>

    <h3>Rows</h3>
    <div class="cnt">
      <div><strong id="vis">0</strong>Shown</div>
      <div><strong id="off">0</strong>Hidden</div>
      <div><strong id="tot">0</strong>Matching</div>
    </div>

    <h3>Display</h3>
    <label class="inl"><input type="checkbox" onchange="cb(this.checked)"> Colourblind palette</label>
    <button class="btn" style="width:100%;margin-top:8px" onclick="reset()">Reset all</button>
  </aside>

  <div class="main">
    <div class="ban i">
      <strong>This file works offline.</strong> Everything it needs is inside it:
      no internet, no server and no DNAInsight installation. It makes no network
      requests of any kind. The card border shows repute, green for generally
      favourable, red for generally unfavourable, grey for unclassified or
      genuinely mixed. Traits and risk scores are always grey on purpose,
      because a trait is not good or bad.
    </div>
    {qc_html}
    <div class="card" style="padding:10px 14px">
      <div id="rep" style="font-size:.86em"></div>
    </div>
    <div id="list"></div>
    <div id="more" style="text-align:center;margin:14px 0"></div>
    {inc_html}
  </div>
</div>

<footer>
  <strong>Not a medical document.</strong> DNAInsight is not a medical device and
  this report is not medical advice. Consumer DNA arrays are not clinical-grade
  tests and cover far less than clinical sequencing, so a negative result here
  does not rule anything out. Do not start, stop or change any medication based
  on this file. Confirm any significant finding with a clinically validated test
  and discuss it with a licensed clinician, pharmacist or genetic counsellor.
  <br><br>
  Magnitude and repute shown here are computed by DNAInsight from CC0 and public
  domain evidence, including CPIC guideline levels, ClinVar review status,
  population frequency and publication counts. They are not the SNPedia values
  of the same name.
  <br><br>
  Generated {_esc(payload["generated_at"])} by DNAInsight
  v{_esc(payload["app_version"])} from {len(payload["findings"])} findings.
</footer>

<script type="application/json" id="payload">{_json_island(payload)}</script>
<script>{_JS}</script>
</body>
</html>
"""
