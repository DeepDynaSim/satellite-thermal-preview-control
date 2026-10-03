"""Verify the delivered manifest and essential study/document counts (standard library only)."""
from pathlib import Path
import hashlib, json, zipfile, xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parent.parent
manifest=ROOT/'SHA256SUMS.txt'
assert manifest.is_file(), 'Run from the extracted delivered package containing SHA256SUMS.txt.'
checked=0
for line in manifest.read_text(encoding='utf-8').splitlines():
    expected,relative=line.split('  ',1)
    target=(ROOT/relative).resolve()
    assert target.is_relative_to(ROOT.resolve()), f'Unexpected manifest path: {relative}'
    assert target.is_file(), f'Missing file: {relative}'
    digest=hashlib.sha256(target.read_bytes()).hexdigest()
    assert digest==expected, f'Checksum mismatch: {relative}'
    checked+=1

data=ROOT/'data'
lhs=json.loads((data/'lhs_report.json').read_text(encoding='utf-8'))
assert lhs['N']==320 and len(lhs['results'])==320
assert all(set(case)=={'LQR','Oracle24','RawNN','FilteredNN'} for case in lhs['results'])
refs=json.loads((data/'references/crossref_verified.json').read_text(encoding='utf-8'))
assert len(refs)==20 and all(r['status']=='verified' for r in refs)
equations=json.loads((data/'equation_inventory.json').read_text(encoding='utf-8'))
assert len(equations)==59
docx=ROOT/'Thermal_Preview_Learning_and_Lyapunov_Safeguarding_Revised.docx'
with zipfile.ZipFile(docx) as archive:
    assert archive.testzip() is None
    xml=ET.fromstring(archive.read('word/document.xml'))
math_count=len(xml.findall('.//{http://schemas.openxmlformats.org/officeDocument/2006/math}oMath'))
assert math_count==144
qa=json.loads((data/'visual_qa.json').read_text(encoding='utf-8'))
assert qa['all_pages_inspected'] and qa['page_count']==30 and qa['defects_remaining']==0
print(json.dumps({'manifest_files_verified':checked,'paired_cases':320,'control_laws_per_case':4,
                  'DOI_records_verified':20,'display_equations':59,'native_Office_Math_objects':math_count,
                  'document_pages':30,'archive_release_status':'See Zenodo_Deposit_Plan.md'},indent=2))
