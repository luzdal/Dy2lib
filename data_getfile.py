# ***************************************************************************************************************************
from dsty_fnc import get_FCF,fnc_idxVm,get_ψNorm,fnc_RvdW,fnc_lam,fnc_slam,val_mass,Evdw_Rvdw,file_c6,fnc_hash,open_file,get_VLJ
# ***************************************************************************************************************************
from dsty_cte import a0,m_e,amu,Eh,h
from dsty_str import STYsí,STYno,STYrarr,STYlarr,STYdarr,STYsp,STYud,STYdt,STYds,STYerr,STYdos,STYsrc,STYset,STYdne
# ***************************************************************************************************************************
import os
path_loc = '.' ; fold_loc = os.listdir(path_loc)
# ***************************************************************************************************************************
import numpy as np
import glob ;
import time

def get_fpath(fpath,fname): return os.path.join(fpath, fname)

def get_fold(path):
    fold = os.listdir(path)
    if fold: return fold
    if not fold: print('-brr')

class MFG_inp():
    def __init__(self,*argv,**kwargs):
        self.argv = argv
        self.frename = {'groundwave.txt':'Wf_g.txt','groundEb.txt':'Eb_g.txt','excitedwave.txt':'Wf_e.txt','excitedEb.txt':'Eb_e.txt'}
        self.inp_lines = {'amu':[1,0],'c6g':[4,0],'c6e':[7,0],'c12g':[4,1],'c12e':[7,1]}

        self.src = 'data_getfile.py | class MFG_inp '

        kv_inp = kwargs     ; self.kv_inp=kv_inp
        kv_dix = {'Ωe':[17,16,15],'iso':164} ;
        kv_dix.update({kk:kwargs[kk] for kk,vv in kv_dix.items() if kk in kv_inp}) ; self.kv_dix=kv_dix
        self.amu = val_mass(**kwargs)
        if kv_dix['iso']==162: self.kk_mm = 'm1m1'
        if kv_dix['iso']==164: self.kk_mm = 'm2m2'

        gΩ,gc12 = 16,361854239 ; self.gΩ,self.gc12 = gΩ,gc12
        kv_gg = {'Ωg':16,'c6g':file_c6(**{'Ωg':16})}
        self.kv_s0 = [{'set':str(eΩ)+'e'+str(gΩ)+'g','c6g':file_c6(**{'Ωg':gΩ}),'c12g':gc12,'c6e':ec6,'c12e':gc12} for eΩ in kv_dix['Ωe'] for ec6 in file_c6(**{'Ωe':eΩ})]

    def file_open(self,*argv,**kwargs):
        print(self.kk_fset)
        g_dix = kwargs['gdix'] ; e_dix = kwargs['edix']
        for kk,vv in self.frename.items():
            f_ind = kwargs['fsearch']+vv
            print(f_ind)
            f_str = f_ind[15:] ; f_pth = get_fpath(path_loc,f_ind)
            str_R = vv+'/R[vdW]' ; str_ψ = vv+'/ψ[vSt]'

    def read_inp(self,*argv,**kwargs):
        '''first call'''
        with open('MFG.inp', 'r') as f: lines=f.readlines()
        now_inp = {}; dix_rowcol = self.inp_lines
        for kk,vv in dix_rowcol.items():
            row,col = vv[0],vv[1]; lin_spl = lines[row].split()[col]; now_inp.update(**{kk:lin_spl})
            if kk in argv: return lin_spl
        if 'dix' in argv: return now_inp

    def write_inp(self,*argv):
        inp_new = self.kv_new; prt_wMFG = '%s>> %s set params' %(STYsp[0:8],inp_new)
        inp_old = self.kv_old; prt_rMFG = '%s%s  %s r/MFG.inp ' %(STYsp[0:8],STYsí,inp_old) if inp_old==inp_new else '%s%s  %s r/MFG.inp' %(STYsp[0:8],STYno,inp_old)
        if inp_old==inp_new:
            if not self.set_exists():
                print('%s %s up-to-date '%(prt_rMFG,STYrarr)) ; print('%s ! run: nohup ./exe > output.txt & '%(prt_wMFG))
                y_n = input(STYsp[0:8]+'>>> ran (n/y)? ')
                if y_n == 'y': self.file_output()

        if inp_old!=inp_new:
            print('%s'%(prt_rMFG)) ; print('%s'%(prt_wMFG))
            '''second call'''
            with open('MFG.inp', 'r') as f: lines=f.readlines()
            for kk,vv_new in inp_new.items():
                vv_old,rc_lin = inp_old[kk],self.inp_lines[kk]
                rc_row,rc_col = rc_lin[0]  ,rc_lin[1]
                c6e,c12e = inp_new['c6e'],inp_new['c12e']
                c6g,c12g = inp_new['c6g'],inp_new['c12g']
                if vv_old != vv_new:
                    kk_str = kk if len(kk)==4 else '.'+kk
                    kk_prt = print('%s >> updating ....%s(%s → %s)'%(STYsp[0:4],kk_str,vv_old,vv_new))
                    vv_now = lines[rc_row].split() ; vv_now[rc_col] = vv_new
                    if kk=='amu': lines[rc_row] = " ".join(vv_now)+"\n"
                    if kk=='c6g' or kk=='c12g': lines[rc_row] = f"{c6g}\t\t{c12g}\n"
                    if kk=='c6e' or kk=='c12e': lines[rc_row] = f"{c6e}\t\t{c12e}\n"
                    with open('MFG.inp', 'w') as f: f.writelines(lines)

    def file_output(self,**kwargs):
        kv_out = self.frename
        for kk_old,kk_txt in kv_out.items():
            kk_new, sí = self.kk_fset+kk_txt , '%s %s |.........'%(STYsp[0:5],STYsí)
            vv_idx, no = 15-len(kk_old)      , '%s %s |.........'%(STYsp[0:5],STYno)
            if kk_old in fold_loc:
                print('%s%s%s: rename %s %s'%(sí,STYdt[:vv_idx],kk_old,STYrarr,kk_new))
                self.file_setname(**{'oldname':kk_old,'newname':kk_new})
            if kk_old not in fold_loc:
                print('%s%s%s: not found'%(no,STYdt[:vv_idx],kk_old))

    def set_exists(self,**kwargs):
        dix = {} ; kv_frename = self.frename
        for idx,kk in enumerate(kv_frename.values()):
            sí = '%s file set %g'%(STYsí,self.idx)
            vv = self.kk_fset+kk
            sp = STYsp[0:35]+kk
            if vv in fold_loc: dix.update({str(idx):vv})
            if vv not in fold_loc: print('%s %s'%(STYerr,vv[:-8])); return

        if len(dix)==4: return dix
    def file_setname(self,**kwargs):
        if 'oldname' and 'newname' in kwargs: os.rename(os.path.join(path_loc, kwargs['oldname']), os.path.join(path_loc, kwargs['newname'])) ; time.sleep(1)

    def run_exists(self,**kwargs):
        run_A = {'set':'runA','c6g':'1.00','c6e':'1.00','c12g':'1.00','c12e':'1.00'}
        run_B = {'set':'runB','c6e':'0.95'} ; run_B.update({kk:vv for kk,vv in run_A.items() if kk not in run_B})
        run_C = {'set':'runC','c6e':'1.05'} ; run_C.update({kk:vv for kk,vv in run_A.items() if kk not in run_C})
        run_D = {'set':'runD','c6e':'0.90'} ; run_D.update({kk:vv for kk,vv in run_A.items() if kk not in run_D})
        run_E = {'set':'runE','c6e':'1.10'} ; run_E.update({kk:vv for kk,vv in run_A.items() if kk not in run_E})
        run_F = {'set':'runF','c6e':'1.10'} ; run_F.update({kk:vv for kk,vv in run_A.items() if kk not in run_F})

        ''' out: list of runs '''
        run_ls = [run_A,run_B,run_C,run_D,run_E,run_F]
        return run_ls
    def run_inp(self,*argv,**kwargs):
        ii=-1
        if self.run_exists():
            arg=argv[0] ; letters=[run['set'][-1] for idx,run in enumerate(self.run_exists())]
            def_src='%s %s | def run_inp'%(STYsrc,self.src)
            if arg not in letters:
                prt_err='%s run%s D.N.E.'%(STYerr,arg)
                prt_dos='%s add new run under <def run_exists> '%(STYdos)  ; print(def_src),print(prt_err),print(prt_dos)#,print()
                return 0

            if arg in letters:
                kk_arg = 'set'
                vv_arg = 'run%s'%(arg)
                kv_run = [kv for kv in self.run_exists() if kv['set'] == vv_arg][0]
                parts = [] ;
                for kk,vv in kv_run.items():
                    add = '(%s|%s); '%(vv,kk)
                    if vv=='1.00': add = '(1|%s); '%(kk)
                    if kk=='set' : add = '%s = (scale|value) = '%(vv)
                    parts.append(add)
                vv_prt='%s %s'%(STYset,''.join(parts))                     ; print(def_src),print(vv_prt)

                kv_inp = {'amu':str(self.amu)+'d0'} ; kv_inp_list,kk_out_list = [],[]
                for vv_idx,kv_s0 in enumerate(self.kv_s0 ):
                    self.idx=vv_idx
                    kv_inp_copy={'Ωg':kv_s0['set'][3:5],'Ωe':kv_s0['set'][0:2]}
                    for kk,sn in kv_run.items():
                        vv_s0 = kv_s0[kk]
                        if kk!='set':
                            sn = float(sn)
                            if 'c6'  in kk:
                                if 'g' in kk: kk_gg = '%.0fg'%vv_s0
                                if 'e' in kk: kk_ee = '%.0fe'%vv_s0
                                kv_inp.update({ kk:'%.6f'%(vv_s0*sn) })
                            if 'c12' in kk: kv_inp.update({ kk:'%sd8'%(vv_s0*sn*1e-8)})
                    kv_inp_copy.update(kv_inp)
                    kv_inp_list.append(kv_inp_copy)

                    kk_fset = '%s[%s_%s]'%(kv_run['set'],self.kk_mm,kk_ee)
                    self.kk_fset = kk_fset

                    check_file = self.set_exists()
                    if check_file:
                        parts = [] ; kk_out = []
                        for kk, vv in check_file.items():
                            kk_out.append(vv)
                            add='%s'%(vv) ; parts.append(add)
                        vv_prt = ','.join(parts) ; #print(sp_prt+'   %s : files = %s'%(vv_idx,vv_prt))
                        kk_out_list.append(kk_out)
                    if not check_file: self.kv_old,self.kv_new=self.read_inp('dix'),kv_inp; self.write_inp(); return

                kv_run.update({kk:float(vv) for kk,vv in kv_run.items() if kk!='set'})
                for idx,kv in enumerate(kv_inp_list):
                    for kk,vv in kv.items():
                        if 'files' in kk: pass
                        if 'd0' in vv: kv.update({kk:float(vv[:-2])})
                        if 'd8' in vv: kv.update({kk:float(vv[:-2])*1e8})
                        if 'd' not in vv: kv.update({kk:float(vv)})
                    kv_out = {'dat':kk_out_list[idx]}
                    kv.update(kv_out)
                return {'run':kv_run,'inp':kv_inp_list}



        #         # return
        #
        # for arg in argv:
        #     ii+=1
        #     if arg == 'A' : kv_run = run_A
        #     if arg == 'B' : kv_run = run_B
        #     if arg == 'C' : kv_run = run_C
        #     if arg == 'D' : kv_run = run_D
        #     if arg == 'E' : kv_run = run_E
        #     else:
        #         print('aaaa')
        #         print('\n \'run\':\'%s\' not yet defined'%arg)
        #
        #     parts = []
        #     for kk, vv in kv_run.items():
        #         add = '(%s|%s); '%(vv,kk)
        #         if vv=='1.00': add = '(1|%s); '%(kk)
        #         if kk=='set' : add = '  %s: run_%s = (scale|value) = '%(kk,vv[-1])
        #         parts.append(add)
        #     sp_prt=STYsp[0:10] ; vv_prt=''.join(parts) ; #print(sp_prt+vv_prt)
        #
        #     kv_inp_list,kv_inp = [],{'amu':str(self.amu)+'d0'}
        #     kk_out_list = []
        #     for vv_idx,kv_s0 in enumerate(self.kv_s0 ):
        #         self.idx=vv_idx
        #         kv_inp_copy={'Ωg':kv_s0['set'][3:5],'Ωe':kv_s0['set'][0:2]}
        #         for kk,sn in kv_run.items():
        #             vv_s0 = kv_s0[kk]
        #             if kk!='set':
        #                 sn = float(sn)
        #                 if 'c6'  in kk:
        #                     if 'g' in kk: kk_gg = '%.0fg'%vv_s0
        #                     if 'e' in kk: kk_ee = '%.0fe'%vv_s0
        #                     kv_inp.update({ kk:'%.6f'%(vv_s0*sn) })
        #                 if 'c12' in kk: kv_inp.update({ kk:'%sd8'%(vv_s0*sn*1e-8)})
        #         kv_inp_copy.update(kv_inp)
        #         kv_inp_list.append(kv_inp_copy)
        #
        #         kk_fset = '%s[%s_%s]'%(kv_run['set'],self.kk_mm,kk_ee)
        #         self.kk_fset = kk_fset
        #
        #         check_file = self.set_exists()
        #         if check_file:
        #             parts = [] ; kk_out = []
        #             for kk, vv in check_file.items():
        #                 kk_out.append(vv)
        #                 add='%s'%(vv) ; parts.append(add)
        #             vv_prt = ','.join(parts) ; #print(sp_prt+'   %s : files = %s'%(vv_idx,vv_prt))
        #             kk_out_list.append(kk_out)
        #         if not check_file: self.kv_old,self.kv_new=self.read_inp('dix'),kv_inp; self.write_inp(); return
        #
        #     kv_run.update({kk:float(vv) for kk,vv in kv_run.items() if kk!='set'})
        #     for idx,kv in enumerate(kv_inp_list):
        #         for kk,vv in kv.items():
        #             if 'files' in kk: pass
        #             if 'd0' in vv: kv.update({kk:float(vv[:-2])})
        #             if 'd8' in vv: kv.update({kk:float(vv[:-2])*1e8})
        #             if 'd' not in vv: kv.update({kk:float(vv)})
        #         kv_out = {'dat':kk_out_list[idx]}
        #         kv.update(kv_out)
        #     return {'run':kv_run,'inp':kv_inp_list}
