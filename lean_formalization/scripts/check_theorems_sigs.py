#!/usr/bin/env python3
# Compare every Lean signature quoted in docs/THEOREMS.md (binders and conclusion, docstrings
# ignored) with the tree. A signature that elides a proof term as `…` or `_` reports as
# NO VERBATIM MATCH only when even the elided form does not match; read it by hand. Usage: python3 scripts/check_theorems_sigs.py docs/THEOREMS.md
import os,re,sys
# the Lean tree, located relative to this script so that the check runs on any machine
ROOT=os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","lean","StackedSVD"))
def norm(s): return re.sub(r'\s+',' ',s).strip()
tree=[]
for dp,dn,fn in os.walk(ROOT):
    for f in fn:
        if f.endswith(".lean"): tree.append((os.path.join(dp,f).replace(ROOT+'/',''),norm(open(os.path.join(dp,f)).read())))
txt=open(sys.argv[1]).read()
blocks=re.findall(r'```lean\n(.*?)```',txt,re.S)
ok=bad=0; seen=set()
for b in blocks:
    # strip docstrings and comments
    b=re.sub(r'/-.*?-/','',b,flags=re.S); b=re.sub(r'--[^\n]*','',b)
    for m in re.finditer(r'(?m)^(?:@\[[^\]]*\]\s*)?(?:noncomputable\s+)?(theorem|def|structure|abbrev|lemma)\s+([^\s(:{\[]+)(.*?)(?=^(?:@\[|noncomputable|theorem|def|structure|abbrev|lemma)\s|\Z)',b,re.S):
        kind,name,body=m.group(1),m.group(2),m.group(3)
        sig=body.split(':=')[0] if kind!='structure' else body.split(' where')[0]
        short=name.split('.')[-1]
        s=norm(f"{kind} {short}{sig}")
        # a declaration written with a dotted name in the tree (`def A.B ...`) is quoted that way
        sfull=norm(f"{kind} {name}{sig}")
        if len(s)<25 or name in seen: continue
        seen.add(name)
        # `…` and a bare `_` stand for an elided proof term; match anything there
        pat=re.escape(s).replace('…','\u2026')
        pat=re.sub(r'\\…|…','.*?',pat)
        pat=re.sub(r'(?<=\\ )_(?=\\ |\\\)|$)','.*?',pat)
        hits=[f for f,t in tree if s in t or sfull in t or re.search(pat,t)]
        if hits: ok+=1
        else:
            bad+=1
            # locate name in tree for diagnosis
            where=[f for f,t in tree if re.search(r'(theorem|def|structure|abbrev) '+re.escape(name)+r'\b',t)]
            tag = "NO VERBATIM MATCH" if where else "NAME NOT IN TREE"
            print(tag+":",name," in:",where[:2] if where else "-")
print(f"signatures {ok+bad}, verbatim {ok}, no verbatim match {bad}")
print("A 'NO VERBATIM MATCH' line means the quoted signature does not appear character for")
print("character in the tree. That is expected for a signature the document abbreviates on")
print("purpose; docs/THEOREMS.md says which ones and why. It is a defect only when the name is")
print("there and the abbreviation is not documented. 'NAME NOT IN TREE' is always a defect.")
