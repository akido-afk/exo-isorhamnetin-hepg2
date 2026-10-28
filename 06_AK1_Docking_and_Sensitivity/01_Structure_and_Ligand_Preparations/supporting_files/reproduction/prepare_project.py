"""Create an independent, portable working copy of the docking project."""
from pathlib import Path
import argparse,json,shutil,os,sys,importlib.util,subprocess
def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--output',type=Path,required=True)
 ap.add_argument('--fresh-docking',action='store_true',help='Keep supplied poses in reference_docking/ and create an empty docking/ directory')
 ap.add_argument('--check-environment',action='store_true',help='Check the current Python environment before creating files')
 args=ap.parse_args();package=Path(__file__).resolve().parents[3];out=args.output.resolve()
 if out.exists():raise SystemExit('Choose a NEW output directory; existing directories are never overwritten.')
 if out==package or out.is_relative_to(package):raise SystemExit('Choose a working directory outside the docking package.')
 if args.check_environment:
  result=subprocess.run([sys.executable,str(Path(__file__).with_name('run_analysis.py')),'--python',sys.executable,'--environment-only'])
  if result.returncode:raise SystemExit('The scientific environment check failed before any project files were created.')
 def long(p):
  s=str(p.resolve());return chr(92)*2+'?'+chr(92)+s if os.name=='nt' and not s.startswith(chr(92)*2+'?'+chr(92)) else s
 def source_file(relative):
  # Resolve the preparation chapter from this script's location.
  relative=relative.replace('01_Structure_and_Ligand_Preparation/',Path(__file__).resolve().parents[2].name+'/',1)
  p=(package/relative).resolve()
  if not p.is_relative_to(package):raise SystemExit('Invalid package file path: '+relative)
  if not Path(long(p)).is_file():raise SystemExit('Required supplied file missing: '+relative)
  return p
 layout=json.loads(Path(__file__).with_name('original_layout.json').read_text(encoding='utf-8'))
 overrides=json.loads(Path(__file__).with_name('portable_sources.json').read_text(encoding='utf-8'))
 # Validate all inputs before creating a partially populated working copy.
 for rel,src in layout.items():
  if not (out/rel).resolve().is_relative_to(out):raise SystemExit('Invalid original-layout path: '+rel)
  source_file(src)
 for name,src in overrides.items():source_file(src)
 vina_source=source_file('01_Structure_and_Ligand_Preparation/supporting_files/software/vina_win.exe')
 out.mkdir(parents=True)
 shutil.copy2(long(Path(__file__).with_name('run_analysis.py')),long(out/'run_analysis.py'))
 vina=out/'software/vina_win.exe';vina.parent.mkdir()
 shutil.copy2(long(vina_source),long(vina))
 for rel,src in layout.items():
  dest=Path(long(out/rel));dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(long(source_file(src)),str(dest))
 for name,src in overrides.items():shutil.copy2(long(source_file(src)),long(out/'structures'/name))
 for source in (out/'structures').iterdir():
  if source.suffix.lower() not in ['.py','.r','.bat','.cmd','.pml']:continue
  code=source.read_text(encoding='utf-8-sig')
  revised=code.replace(r'C:\Users\akidi\Desktop\ISO_AK_Docking',out.as_posix()).replace('C:/Users/akidi/Desktop/ISO_AK_Docking',out.as_posix())
  revised=revised.replace(r'C:\Users\akidi\Desktop\Docking tools\vina_win.exe',vina.as_posix())
  if revised!=code:source.write_text(revised,encoding='utf-8')
 if args.fresh_docking:
  if (out/'reference_docking').exists():raise SystemExit('Unexpected reference_docking directory')
  (out/'docking').rename(out/'reference_docking');(out/'docking').mkdir()
 print('Working copy prepared:',out)
 print('Run scripts from:',out/'structures')
 print('Python interpreter:',sys.executable)
 print('No calculations were launched.')
 if args.fresh_docking:print('New search outputs: docking/. Supplied reference poses: reference_docking/.')
 else:print('Saved docking is retained. Completed outputs are skipped by the portable batches; use --fresh-docking for a new search.')
if __name__=='__main__':main()
