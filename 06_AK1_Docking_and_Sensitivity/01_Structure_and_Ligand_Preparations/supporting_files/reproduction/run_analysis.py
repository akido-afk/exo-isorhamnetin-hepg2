"""Run an original docking-analysis script with its working directory and runtime configured."""
from pathlib import Path
import argparse,os,sys,subprocess,shutil,tempfile,re

def scientific_env(python):
 prefix=Path(python).resolve().parent;env=dict(os.environ)
 if (prefix/'Library/bin').is_dir():
  dirs=[prefix,prefix/'Library/mingw-w64/bin',prefix/'Library/usr/bin',prefix/'Library/bin',prefix/'Scripts',prefix/'bin']
  env['PATH']=os.pathsep.join(str(p) for p in dirs if p.exists())+os.pathsep+env.get('PATH','')
  env['CONDA_PREFIX']=str(prefix);env['CONDA_DEFAULT_ENV']=prefix.name
  env['CONDA_DLL_SEARCH_MODIFICATION_ENABLE']='1'
 env['PYTHONIOENCODING']='utf-8';env['PYTHONDONTWRITEBYTECODE']='1'
 return env

def find_program(explicit,candidates,label):
 if explicit:
  p=Path(explicit)
  if p.is_file():return str(p.resolve())
  q=shutil.which(explicit)
  if q:return q
  raise SystemExit(label+' executable not found: '+explicit)
 for p in candidates:
  if p and Path(p).is_file():return str(Path(p).resolve())
 raise SystemExit(label+' is not available. Install it and specify its executable path using --'+label.lower()+'.')

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('script',nargs='?');ap.add_argument('--project',type=Path,default=Path(__file__).resolve().parent)
 ap.add_argument('--python',default=sys.executable);ap.add_argument('--rscript');ap.add_argument('--pymol');ap.add_argument('--r-library',type=Path)
 ap.add_argument('--allow-docking',action='store_true');ap.add_argument('--environment-only',action='store_true');ap.add_argument('--preview',action='store_true',help='Render a separate PyMOL preview at up to 640 pixels')
 args=ap.parse_args();python=find_program(args.python,[],'Python');env=scientific_env(python)
 if args.environment_only:
  code="import importlib.util,sys;mods=['numpy','pandas','scipy','rdkit','xlsxwriter','meeko','gemmi','pdbfixer'];missing=[m for m in mods if importlib.util.find_spec(m) is None];print('Python:',sys.executable,flush=True);print('Missing packages:',', '.join(missing) or 'none',flush=True);sys.exit(2) if missing else None;import numpy as n;n.linalg.svd(n.eye(3));print('Numerical runtime ready',flush=True)"
  r=subprocess.run([python,'-X','utf8','-c',code],env=env)
  if r.returncode:raise SystemExit('The selected scientific environment is incomplete or cannot load its numerical libraries. Use the docking environment and install its declared dependencies with the same Python.')
  return
 if not args.script:ap.error('Specify an original source filename')
 project=args.project.resolve();structures=project/'structures';script=structures/args.script
 if not script.resolve().is_relative_to(structures) or not script.is_file():raise SystemExit('Script not found in the prepared project structures/: '+args.script)
 suffix=script.suffix.lower()
 if suffix in ['.bat','.cmd'] or script.name=='run_AK1_multiseed.py':
  if not args.allow_docking:raise SystemExit('This is a docking launcher. Use --allow-docking only when you intend to run it. For new searches first prepare a --fresh-docking working copy.')
 temporary=None;expected_images=[]
 if suffix=='.py':command=[python,'-X','utf8',str(script)]
 elif suffix=='.r':
  candidates=[shutil.which('Rscript')]
  rhome=env.get('R_HOME')
  if rhome:candidates.extend([Path(rhome)/'bin/Rscript.exe',Path(rhome)/'bin/x64/Rscript.exe'])
  candidates.extend(sorted(Path(env.get('ProgramFiles',r'C:\Program Files')).glob('R/R-*/bin/Rscript.exe'),reverse=True))
  rscript=find_program(args.rscript,candidates,'Rscript')
  # Use the selected R version's user library; do not install packages implicitly.
  candidates=[Path.home()/'AppData/Local/R/win-library']
  if args.r_library:env['AK1_R_LIBRARY']=str(args.r_library.resolve())
  env['AK1_R_USER_BASE']=str(candidates[0])
  for key in ['LC_ALL','LC_CTYPE','LC_COLLATE','LC_MONETARY','LC_TIME','LANG']:
   if env.get(key) in ['C.UTF8','C.UTF-8']:env.pop(key)
  bridge="args<-commandArgs(trailingOnly=TRUE);if(.Platform$OS.type=='windows')suppressWarnings(try(Sys.setlocale('LC_CTYPE','English_United States.utf8'),silent=TRUE));options(encoding='UTF-8');ver<-paste(R.version$major,strsplit(R.version$minor,'.',fixed=TRUE)[[1]][1],sep='.');libs<-c(Sys.getenv('AK1_R_LIBRARY'),file.path(Sys.getenv('AK1_R_USER_BASE'),ver),file.path(Sys.getenv('LOCALAPPDATA'),'R','win-library',ver));libs<-libs[nzchar(libs)&dir.exists(libs)];.libPaths(unique(c(libs,.libPaths())));install.packages<-function(...)stop('Install the required R packages explicitly before running this source');write.csv<-function(x,file='',...,fileEncoding='UTF-8')utils::write.csv(x=x,file=file,...,fileEncoding=fileEncoding);ggsave<-function(filename,...,device=NULL){if(is.null(device)&&grepl('[.]pdf$',filename,ignore.case=TRUE)&&capabilities('cairo'))device<-grDevices::cairo_pdf;ggplot2::ggsave(filename=filename,...,device=device)};source(args[1],chdir=FALSE,encoding='UTF-8')"
  command=[rscript,'--vanilla','-e',bridge,str(script)]
 elif suffix=='.pml':
  prefix=Path(python).parent
  pymol=find_program(args.pymol,[shutil.which('pymol'),prefix/'Scripts/pymol.exe',prefix/'bin/pymol'],'PyMOL')
  text=script.read_text(encoding='utf-8-sig')
  if args.preview:
   lines=[]
   for line in text.splitlines():
    if line.strip().startswith('mset '):break
    size=re.match(r'^(\s*)(ray|viewport)\s+(\d+)\s*,\s*(\d+)',line)
    if size:
     w,h=int(size[3]),int(size[4]);scale=min(1.0,640/max(w,h));line=f'{size[1]}{size[2]} {max(1,round(w*scale))}, {max(1,round(h*scale))}'
    if line.strip().startswith('png '):
     raw=line.strip()[4:].split(',')[0].strip().strip('"');dest=project/'preview_figures'/script.stem/Path(raw.replace('\\','/')).name
     line='png '+dest.as_posix()+(','+line.split(',',1)[1] if ',' in line else '')
    lines.append(line)
   text='\n'.join(lines)+'\nquit\n'
   with tempfile.NamedTemporaryFile(mode='w',encoding='utf-8',suffix='.pml',dir=structures,delete=False) as f:f.write(text);temporary=Path(f.name)
  for line in text.splitlines():
   line=line.strip()
   if line.startswith('load '):
    raw=line[5:].split(',')[0].strip().strip('"');p=structures/Path(raw.replace('\\','/'))
    if not p.is_file():raise SystemExit('Required PyMOL input missing: '+str(p))
   if line.startswith(('png ','mpng ')):
    raw=line.split(' ',1)[1].split(',')[0].strip().strip('"');dest=structures/Path(raw.replace('\\','/'));dest.parent.mkdir(parents=True,exist_ok=True)
    if line.startswith('png '):expected_images.append((dest,dest.stat().st_mtime_ns if dest.exists() else None))
  command=[pymol,'-cq',str(temporary or script)]
 elif suffix in ['.bat','.cmd']:command=[env.get('COMSPEC','cmd.exe'),'/d','/c',str(script)]
 else:raise SystemExit('Supported sources: .py, .R, .pml, .bat and .cmd')
 print('Running:',script.name,flush=True)
 if suffix in ['.bat','.cmd']:
  r=subprocess.run(command,cwd=structures,env=env,input='\n',text=True)
 else:r=subprocess.run(command,cwd=structures,env=env)
 if temporary and temporary.exists():temporary.unlink()
 if r.returncode:
  if r.returncode in [3228369023,-1066598273]:raise SystemExit('The numerical runtime could not load required DLLs. Run --environment-only with the selected scientific Python to check the environment.')
  raise SystemExit(r.returncode if 0<r.returncode<256 else 1)
 for image,previous_mtime in expected_images:
  if not image.is_file() or image.stat().st_size<24 or image.stat().st_mtime_ns==previous_mtime:raise SystemExit('PyMOL did not regenerate expected image: '+str(image))
  with image.open('rb') as f:
   if f.read(8)!=b'\x89PNG\r\n\x1a\n':raise SystemExit('Invalid PNG output: '+str(image))
if __name__=='__main__':main()
