from data_getfile import MFG_inp
# ***************************************************************************************************************************
from dsty_fnc import file_path,file_rename,get_FCF,fnc_idxVm,get_ψNorm,fnc_RvdW,fnc_lam,fnc_slam,val_mass,Evdw_Rvdw,file_c6,fnc_hash,open_file,get_VLJ
from dsty_str import STYsí, STYno, STYrarr, STYlarr, STYdarr, STYsp, STYud, STYdt, STYds, STYerr, STYdos, STYsrc,STYinp
# ***************************************************************************************************************************
import os
path_loc = '.' ; fold_loc = os.listdir(path_loc)
# ***************************************************************************************************************************
import numpy as np
import datetime

class MFG_out():
    def __init__(self,*argv,**kwargs):

        self.kk_f = 'fname'
        self.iso  = kwargs['iso']

        self.src = 'data_getvals.py | class MFG_out '#%(STYrarr,STYrarr)

        data = MFG_inp(**{'iso':self.iso}).run_inp(kwargs['run'])
        self.kv_init = None

        def_src='%s %s | def __init__ '%(STYsrc,self.src)
        if data==0:
            prt_err = '%s output from  <MFG_inp.run_inp> = 0'%(STYerr)     ; print(def_src,'\n'+prt_err)

        if data!=0:
            if data == None:
                # pass
                prt_err = '%s no files found for dataset ''run%s'''%(STYerr,kwargs['run'])

                # print(def_src,'\n'+prt_err)

            else:

                data_files = [kv['dat'] for idx,kv in enumerate(data['inp'])]
                data_inp = ' | '.join(data_files[0])                                ; print(def_src), print('%s %s'%(STYinp,data_inp))

                self.data    = data
                self.run_kv  = {'s%s'%(kk[1:]):vv for kk,vv in self.data ['run'].items()}
                self.data_set = self.data ['inp']
                self.name_set = [kv['dat'][0][:-8] for idx,kv in enumerate(self.data_set)]
                self.data_len = len(self.data_set)
                self.data_ran = range(self.data_len)

                amu = [kv['amu'] for idx,kv in enumerate(self.data_set)] ; self.amu = sum(amu)/len(amu)

                kv_prt = {} ; kv_prt.update({'arg':arg for arg in argv if arg=='prt'})
                kv_dim = {} ; kv_dim.update({'dim':arg for arg in argv if arg=='dim'})

                if kv_prt: self.fdata_kv(kv_prt['arg'])
                else:
                    self.data_kv = self.fdata_kv()
                    self.edim = None
                    if kv_dim:
                        self.kv_init = self.fdata_kv(kv_dim['dim'])
                        # self.edim    = self.fdata_kv(kv_prt['dim'])
                    call = self.__str__()
                    # self.edim = self.fdata_kv(kv_prt['dim']) if kv_dim else None

    def __str__(self,*argv):
        fdata_kv = self.fdata_kv()

        if fdata_kv:
            for idx, kv in enumerate(fdata_kv):
                g6 ,e6 = '%.2f'%kv['c6g'], '%.2f'%kv['c6e']; s6 ='%sc6%s'%(STYud[0:3],STYud[0:4])
                g12,e12='%.4e'%kv['c12g'],'%.4e'%kv['c12e']; s12='%sc12%s'%(STYud[0:4],STYud[0:5])

                gψ='(%g)'%kv['gψdat'][1].shape[0]; gV0='(%.3e)'%min(kv['gV[Eh]']); gE0='(%.3e)'%min(kv['gE[Eh]'])
                eψ='(%g)'%kv['eψdat'][1].shape[0]; eV0='(%.3e)'%min(kv['eV[Eh]']); eE0='(%.3e)'%min(kv['eE[Eh]'])
                sψ ='%sdimψ%s'%(STYud[0:1],STYud[0:1])
                sE0='%sminEb[Eh]%s'%(STYud[0:2],STYud[0:3])
                sV0='%sminV[Eh]%s'%(STYud[0:3],STYud[0:3])
                gR,eR= '%.3e'%kv['gRvdw']  ,'%.3e'%kv['eRvdw']  ; sR='%svdwR%s'%(STYud[0:3],STYud[0:4])
                gE,eE= '%.3e'%kv['gEvdw']  ,'%.3e'%kv['eEvdw']  ; sE='%svdwE%s'%(STYud[0:3],STYud[0:4])
                fiso = 'iso%s'%(self.iso)#kv['run'][-4:]
                frun = 'run%s'%(kv['run'][3])

                t_prt = '~%s/%s'%(fiso,STYds[0:83])
                s_prt = '~%s|%s.%s.%s.%s.%s.%s.%s|~%s'%(frun,sψ,sE0,sV0,s6,s12,sR,sE,fiso)
                g_prt = '%s~g| %s | %s | %s | %s | %s | %s | %s |'%(STYsp[0:3],gψ,gE0,gV0,g6,g12,gR,gE)
                e_prt = '%s~e| %s | %s | %s | %s | %s | %s | %s |'%(STYsp[0:3],eψ,eE0,eV0,e6,e12,eR,eE)

                if idx==0: print(s_prt,'\n'+g_prt)
                print(e_prt)

    def fdata_kv(self,*argv):
        kv_fd_list = []
        for idx in self.data_ran:
            kv_data,self.run_kv['set'] = self.data_set[idx].copy(),self.name_set[idx]
            kv_data.pop('dat',None) ; kv_data.update({'idx':idx}) ;
            kv_g = {'c6g':kv_data['c6g'],'amu':kv_data['amu']}; kv_data.update(Evdw_Rvdw(**kv_g))
            kv_e = {'c6e':kv_data['c6e'],'amu':kv_data['amu']}; kv_data.update(Evdw_Rvdw(**kv_e))
            kv_data.update(self.run_kv)
            kv_data.update({'run':kv_data['set'][0:4]+'/'+kv_data['set'][5:9]})

            g_Evdw = [vv for kk,vv in Evdw_Rvdw(**kv_g).items() if 'gE' in kk][0]
            e_Evdw = [vv for kk,vv in Evdw_Rvdw(**kv_e).items() if 'eE' in kk][0]
            g_Rvdw = [vv for kk,vv in Evdw_Rvdw(**kv_g).items() if 'gR' in kk][0]
            e_Rvdw = [vv for kk,vv in Evdw_Rvdw(**kv_e).items() if 'eR' in kk][0]

            kv_Wf = self.data_Wf()[idx].copy()
            kv_Wf.update({'gR[a0]':vv*g_Rvdw for kk,vv in kv_Wf.items() if 'gR' in kk}); kv_Wf.pop('gRdat',None)
            kv_Wf.update({'eR[a0]':vv*e_Rvdw for kk,vv in kv_Wf.items() if 'eR' in kk}); kv_Wf.pop('eRdat',None)
            kv_Eb = self.data_Eb()[idx].copy()
            kv_Eb.update({'gE[Eh]':vv*g_Evdw for kk,vv in kv_Eb.items() if 'gE' in kk}); kv_Eb.pop('gEdat',None)
            kv_Eb.update({'eE[Eh]':vv*e_Evdw for kk,vv in kv_Eb.items() if 'eE' in kk}); kv_Eb.pop('eEdat',None)
            kv_data.update(**kv_Wf),kv_data.update(**kv_Eb)

            kv_data.update(self.VLJ_kv(**kv_data))#; kv_data.update({'eψdim'})
            kv_data.update(self.FCF_kv(**kv_data))#; kv_data.update({'eψdim'})
            # [print(kk) for kk,vv in kv_data.items() if idx==0]
            kv_fd_list.append(kv_data)

        if kv_fd_list:
            if 'dim' in argv:
                dim_e = [kv['eψdat'][1].shape[0]  for idx,kv in enumerate(kv_fd_list)]
                dim_g = [kv['gψdat'][1].shape[0]  for idx,kv in enumerate(kv_fd_list)]
                return {'gdim':dim_g,'edim':dim_e}
            else : return kv_fd_list


    def FCF_kv(self,*argv,**kwargs):
        gψnor = get_ψNorm('g',**{'rpts':kwargs['gR[a0]'],'vpts':kwargs['gψdat']})
        eψnor = get_ψNorm('e',**{'rpts':kwargs['eR[a0]'],'vpts':kwargs['eψdat'],'rpts_g':kwargs['gR[a0]']})
        return {'FCFmat':get_FCF(**{'rpts_g':kwargs['gR[a0]'],'ψNor_g':gψnor,'ψNor_e':eψnor})}

    def VLJ_kv(self,*argv,**kwargs):
        def fnc_VLJ(**kv):
            dix = {xx:vv for kk,vv in kv.items() for xx in ['c6','c12','R'] if xx in kk}
            return  (dix['c12'] / (dix['R']**12)) - (dix['c6'] / (dix['R']**6))
        kv_ge   = {kk:vv for kk,vv in kwargs.items() for xx in ['c6','c12','R'] if xx in kk}
        gkv,ekv = {kk:vv for kk,vv in kv_ge.items() if 'g' in kk},{kk:vv for kk,vv in kv_ge.items() if 'e' in kk}
        kv      = {'gV[Eh]':fnc_VLJ(**gkv),'eV[Eh]':fnc_VLJ(**ekv)}
        return kv

    def data_Wf(self,*argv):
        kv_list = []
        for idx,kv_dat in enumerate(self.data_set):
            fname_gWf = [f for f in kv_dat['dat'] if 'Wf' in f if 'g.txt' in f][0]
            fname_eWf = [f for f in kv_dat['dat'] if 'Wf' in f if 'e.txt' in f][0]
            kv_Rψ  = {'set':kv_dat['dat'][0][:-8]}
            kv_eWf = self.readf_Wf(*argv,**{'fname':fname_eWf}) ; kv_Rψ.update(kv_eWf)
            kv_gWf = self.readf_Wf(*argv,**{'fname':fname_gWf}) ; kv_Rψ.update(kv_gWf)
            kv_list.append(kv_Rψ)
        if kv_list:
            parts = ['%s,'%kk for kk,vv in kv_list[0].items()]
            if argv: print(' %s data_Wf() | [k_o=%s]'%(STYrarr,''.join(parts)[:-1]))
            return kv_list

    def data_Eb(self,*argv):
        kv_list = []
        for idx,kv_dat in enumerate(self.data_set):
            fname_gEb = [f for f in kv_dat['dat'] if 'Eb' in f if 'g.txt' in f][0]
            fname_eEb = [f for f in kv_dat['dat'] if 'Eb' in f if 'e.txt' in f][0]
            kv_VE  = {'set':kv_dat['dat'][0][:-8]}
            kv_eEb = self.read_Eb(*argv,**{'fname':fname_eEb}) ; kv_VE.update(kv_eEb)
            kv_gEb = self.read_Eb(*argv,**{'fname':fname_gEb}) ; kv_VE.update(kv_gEb)
            kv_list.append(kv_VE)
        if kv_list:
            parts = ['%s,'%kk for kk,vv in kv_list[0].items()]
            if argv: print(' %s data_Eb() | [k_o=%s]'%(STYrarr,''.join(parts)[:-1]))
            return kv_list

    def read_Eb(self,*argv,**kwargs):
        kk_i = '   %s read_Eb(**{k_i=%s:v_i=%s})'%(STYrarr,self.kk_f,kwargs[self.kk_f])
        kv = {}
        for kk,vv in kwargs.items():
            if 'Eb' in vv:
                f_name = vv ; Eb_eg = f_name[-5:-4]
                f_path = file_path(**{'filepath':path_loc,'filename':f_name})
                if os.path.isfile(f_path):
                    Eb_data = np.loadtxt(f_path,skiprows=1)
                    Epts = Eb_data[:,0]
                    kv.update({Eb_eg+'Edat':Epts})
        if kv:
            kk_o = ['%s,'%kk for kk,vv in kv.items()]
            if argv: print('%s | [k_o=%s]'%(kk_i,''.join(kk_o)[:-1]))
            return kv

    def readf_Wf(self,*argv,**kwargs):
        kk_i = '   %s read_Wf(**{k_i=%s:v_i=%s})'%(STYrarr,self.kk_f,kwargs[self.kk_f])
        kv = {}
        for kk,vv in kwargs.items():
            if 'Wf' in vv:
                f_name = vv ; Wf_eg = f_name[-5:-4]
                f_path = file_path(**{'filepath':path_loc,'filename':f_name})
                if os.path.isfile(f_path):
                    Wf_data = np.loadtxt(f_path) ; Rpts,ψpts = Wf_data[:,0],Wf_data[:,1:]
                    kv.update({Wf_eg+'Rdat':Rpts,Wf_eg+'ψdat':ψpts})
        if kv:
            kk_o = ['%s,'%kk for kk,vv in kv.items()]
            if argv: print('%s | [k_o=%s]'%(kk_i,''.join(kk_o)[:-1]))
            return kv








# for run in ('A'):
#     for iso in [162]:
# # for run in ('A','B','C','D','E'):
# #     for iso in [162,164]:
#         print()
#         file = MFG_out(**{'iso':iso, 'run':run}).__str__()
        # file_kv = file.fdata_kv()
        # if file_kv:
        #     print(len(file_kv))





                # file_kv = file.fdata_kv()
                # if file_kv:
                #     print(len(file_kv))

        #
        # class data_analys():
        #     def __init__(self,*argv,kv_plt,kv_inp):
        #
        #         self.kv_plt = kv_plt
        #         print('\n ↓|| data_call.py || ')
        #
        #         k_arg = ['plt','pdf','png']
        #         k_plt = argv[0] if argv and argv[0] in k_arg else ()
        #
        #         m_arr = [kv_inp['iso']] if 'iso' in kv_inp else [162,164]
        #         s_arr = [kv_inp['s12']] if 's12' in kv_inp else ['s0', 's1', 's2']
        #
        #         for ss in s_arr:
        #             ''' % ---------- % '''
        #             self.mlen,self.dlen = len(m_arr),6 ; ax_g,ax_e,fig = (), (), ()
        #
        #             self.kv_fig = self.fig_dim() ; kv_fig = self.kv_fig
        #             if k_plt:
        #                 fig=plt.figure(figsize=(kv_fig['figW'], kv_fig['figH']), dpi=100)
        #                 gs=gridspec.GridSpec(kv_fig['rowNr'],kv_fig['colNr']+1,figure=fig,height_ratios=kv_fig['rtoH'],width_ratios=kv_fig['rtoW'])
        #                 fig.subplots_adjust(hspace=kv_fig['sepY'],wspace=kv_fig['sepX'])
        #                 self.fig = fig
        #
        #             for midx,mm in enumerate(m_arr):
        #                 d_arr = file_MFG_inp(ss,**{'iso':mm, 's12':ss}).set_dix
        #                 if d_arr:
        #                     for didx,dd in enumerate(d_arr):
        #                         kv_dd={'midx':midx,'didx':didx} ; kv_dd.update(**dd,**kv_plt) ; self.kv_dd=kv_dd
        #
        #                         gΩ,gc6,gc12 = kv_dd['gΩ'],kv_dd['gc6'],kv_dd['gc12']
        #                         eΩ,ec6,ec12 = kv_dd['eΩ'],kv_dd['ec6'],kv_dd['ec12']
        #                         # self.__str__()
        #                         lab_g = r'$%s(%s=%.0f)%s_{%s %s=%g}$'%(mth_Vlj,mth_c6,gc6,mth_vB,mth_ss,mth_gΩ,gΩ)
        #                         lab_e = r'$%s(%s=%.0f)%s_{%s %s=%g}$'%(mth_Vlj,mth_c6,ec6,mth_vB,mth_ss,mth_eΩ,eΩ)
        #
        #                         ttl = kv_dd['ttl']  ;
        #                         sup = lab_g         ; self.labTex = {'cbar':tex_logFCF,'ψg':tex_vg,'ve':tex_ve,'ttl':ttl,'sup':sup,'leg':lab_e}
        #                         ''' % ---------- % '''
        #                         labTex=self.labTex
        #                         if fig:
        #
        #                             if didx==0:
        #                                 ax_e=fig.add_subplot(gs[didx,midx])
        #                                 ''' % subplot title % '''
        #                                 ax_e.set_title(labTex['ttl'],x=kv_fig['ttlX'],y=kv_fig['ttlY'])
        #                             else : ax_e=fig.add_subplot(gs[didx,midx],sharex=ax_e,sharey=ax_e)
        #                             if didx!=self.dlen-1: ax_e.tick_params(axis='x',labelbottom=False,bottom=True)
        #
        #                             if 'fcf' in kv_plt:
        #                                 self.vv_fcf = kv_plt['fcf']
        #                                 self.fcf_call(ax=ax_e)
        #
        #                 ''' % ---------- % '''
        #         # if fig:
        #         #     fig.suptitle(labTex['sup'],x=kv_fig['supX'],y=kv_fig['supY'],fontsize=11,fontweight='bold',color=clr_g)
        #         #     # fig.suptitle(labTex['sup'],x=kv_fig['supX'],y=kv_fig['supX'],fontsize=12.5,fontweight='bold')
        #         #     ''' % ---------- % '''
        #         #     if k_plt and k_plt=='plt': plt.show()
        #         #     if k_plt and k_plt=='png' or k_plt=='png':
        #         #         plt_name = str_dati+str(kv_plt['ve']) ; plt_fold = 'data_runB' ; outp_dir = plt_fold
        #         #         os.makedirs(outp_dir, exist_ok=True)  # creates it if needed
        #         #         plt.savefig(os.path.join(outp_dir,plt_name),bbox_inches='tight',facecolor='white',edgecolor='white')
        #         #         # plt.savefig(os.path.join(outp_dir,str_set),bbox_inches='tight',transparent=True)
        #
        #     def __str__(self,*argv):
        #         [print(kk) for kk,vv in self.kv_dd.items() if self.kv_dd['midx'] and self.kv_dd['didx']==0]
        #
        #
        #     def fcf_call(self,*argv,ax):
        #         labTex = self.labTex
        #         midx = self.kv_dd['midx'] ; kv_dd = self.kv_dd
        #         didx = self.kv_dd['didx'] ; kv_plt = self.kv_plt
        #         fig  = self.fig           ; kv_fig = self.kv_fig
        #
        #         evHi,evLo,gvHi,gvLo = (),(),(),()
        #         lab_x,lab_y = labTex['ψg'], labTex['ve']
        #         dim_g,dim_e = len(kv_dd['FCF'][0,:]),len(kv_dd['FCF'][:,0])
        #
        #         if not isinstance(kv_plt['ve'], str) and len(kv_plt['ve'])==2:evLo,evHi=kv_plt['ve'][0],kv_plt['ve'][-1]
        #         elif isinstance(kv_plt['ve'], str): evLo,evHi=-dim_e,-1
        #         if not isinstance(kv_plt['ψg'], str) and len(kv_plt['ψg'])==2:gvLo,gvHi=kv_plt['ψg'][0],kv_plt['ψg'][-1]
        #         elif isinstance(kv_plt['ψg'], str): gvLo,gvHi=-dim_g,-1
        #
        #         if not evHi and not gvHi  : mat_FCF=kv_dd['FCF']
        #         elif evHi and not gvHi    : mat_FCF=kv_dd['FCF'][evLo:,:] if evHi ==-1 else kv_dd['FCF'][evLo:evHi+1,:]
        #         elif not evHi and gvHi    : mat_FCF=kv_dd['FCF'][:,gvLo:] if gvHi ==-1 else kv_dd['FCF'][:,gvLo:gvHi+1]
        #         elif evHi and gvHi        :
        #             if evHi == -1 and gvHi == -1: mat_FCF=kv_dd['FCF'][evLo:,gvLo:]
        #             elif evHi == -1              : mat_FCF=kv_dd['FCF'][evLo:,gvLo:gvHi+1]
        #             elif gvHi == -1              : mat_FCF=kv_dd['FCF'][evLo:evHi+1,gvLo:]
        #             else                          : mat_FCF=kv_dd['FCF'][evLo:evHi+1,gvLo:gvHi+1]
        #
        #         vgLen = mat_FCF.shape[1]; vgRan = np.arange(vgLen); self.vgLen,self.vgRan = vgLen,vgRan
        #         veLen = mat_FCF.shape[0]; veRan = np.arange(vgLen); self.veLen,self.veRan = veLen,veRan
        #
        #
        #         if vgLen<=2 : tix_g = vgRan
        #         else        : tix_g = [x for x in [vgRan[1],vgRan[-2]] ]
        #         if 'tixg' in kv_plt: tix_g = kv_plt['tixg']
        #
        #         self.FCF_mtrx = mat_FCF
        #
        #         self.gvLo,self.gvHi=gvLo,gvHi
        #         self.evLo,self.evHi=evLo,evHi
        #
        #         if self.vv_fcf=='plt': self.fcf_plt(ax=ax)
        #         if self.vv_fcf=='img': self.fcf_img(ax=ax)
        #
        #     def fcf_plt(self,*argv,ax):
        #         labTex = self.labTex
        #         kv_fig= self.kv_fig
        #         kv_dd=self.kv_dd ; midx=kv_dd['midx']; didx=kv_dd['didx']
        #         FCF_mtrx = self.FCF_mtrx
        #         FCF_vals = np.concatenate([mat.ravel() for mat in FCF_mtrx])
        #         FCF_pstv = FCF_vals[FCF_vals > 0]
        #         ψmin, ψmax = max(1e-6, FCF_pstv.min()), FCF_pstv.max()
        #
        #         # Log-transform the FCF matrix
        #         # mat_FCF = np.log(np.clip(kv_dd['FCF'], ψmin, ψmax))  # Ensure no log(0)
        #         mat_FCF = np.log(np.clip(self.FCF_mtrx, ψmin, ψmax))  # Correct: use zoomed-in matrix
        #
        #         gvLo,gvHi=self.gvLo,self.gvHi
        #         evLo,evHi=self.gvLo,self.gvHi
        #
        #         len_g = mat_FCF.shape[1]
        #         len_e = mat_FCF.shape[0]
        #         ran_g = np.arange(len_g)
        #         ran_e = np.arange(len_e)
        #
        #         vg_range = np.arange(mat_FCF.shape[1])  # x-axis: ground state v''
        #         fcf_vals = mat_FCF[0]                   # single row (only 1 excited state ve)
        #
        #         for xidx, ypts in enumerate(fcf_vals):
        #             ax.scatter(xidx,ypts,s=mSz-1,edgecolors='none',color=c_sns[xidx%len(c_sns)])
        #         # labb = f"log({self.labTex['ve']})"
        #         labb = r"$\log\left(|\langle \psi_{v'%g} | \psi_{v''} \rangle|^2\right)$"%(evHi)
        #         ax.plot(vg_range, fcf_vals, '--', color='gray', linewidth=lWd,label=labb)
        #
        #         tix_y_vals = [1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1]
        #         tix_y = np.log(tix_y_vals)
        #
        #         ax.set_yticks(tix_y)
        #         ax.set_yticklabels([f"$10^{{{int(np.log10(y))}}}$" for y in tix_y_vals])
        #
        #         # Set global y-limits
        #         ystp = .75
        #         ax.set_ylim([-14-ystp, ystp])
        #
        #         if midx==0:
        #             ax.tick_params(axis='y',direction='in',left=False,labelleft=False,labelright=False,right=True)
        #             ax.yaxis.set_label_position('right')
        #             ax_pat = Patch(facecolor='none',edgecolor='none',label=labTex['leg'])
        #             ax.legend(handles=[ax_pat],loc='upper left',bbox_to_anchor=(kv_fig['patX'],kv_fig['patY']),frameon=True, fancybox=False,labelcolor=kv_fig['patC'],framealpha=0.0)
        #
        #
        #         else:
        #             ax.legend()
        #             ax.set_yticklabels([f"$10^{{{int(np.log10(y))}}}$" for y in tix_y_vals])
        #             # ax.set_ylabel(f"log({self.labTex['ve']})", color='black')
        #             ax.tick_params(axis='y',direction='in',left=False,labelleft=False,labelright=True,right=True)
        #             ax.yaxis.set_label_position('right')
        #
        #         # tix_x = [x+1 for x in ran_g if x % 20 == 0]
        #         tix_x = [x for x in [ran_g[1],ran_g[-2]] ]
        #
        #         ax.set_xticks(tix_x)
        #         if didx==self.dlen-1:
        #             ax.set_xlabel(self.labTex['ψg'], color='black')
        #             # ax.set_xticklabels([f'{-1-x}' for x in ax.get_xticks()])
        #             ax.set_xticklabels([f'{x-abs(gvLo)}' for x in ax.get_xticks() ])
        #         # if didx!=self.dlen-1: ax.tick_params(axis='x',direction='in',top=False,labeltop=False,labelbottom=False,bottom=False)
        #
        #     def fcf_img(self,*argv,ax):
        #         FCF_mtrx  = self.FCF_mtrx
        #         FCF_vals = np.concatenate([mat.ravel() for mat in FCF_mtrx])
        #         FCF_pstv = FCF_vals[FCF_vals > 0]
        #         ψmin,ψmax= max(1e-6,FCF_pstv.min()),FCF_pstv.max()
        #         norm = colors.LogNorm(vmin=ψmin, vmax=ψmax)
        #         ########################
        #         labTex = self.labTex
        #         midx = self.kv_dd['midx'] ; kv_dd = self.kv_dd
        #         didx = self.kv_dd['didx'] ; kv_plt = self.kv_plt
        #         fig  = self.fig           ; kv_fig = self.kv_fig
        #
        #         lab_x,lab_y = labTex['ψg'], labTex['ve']
        #
        #         mat_FCF=self.FCF_mtrx
        #         gvLo,gvHi=self.gvLo,self.gvHi
        #         evLo,evHi=self.evLo,self.evHi
        #
        #         len_g = len(mat_FCF[0,:]) ; ran_g = np.arange(0,len_g)
        #         len_e = len(mat_FCF[:,0]) ; ran_e = np.arange(0,len_e)
        #
        #         if len_g<=2 : tix_x = ran_g
        #         else        : tix_x = [x for x in [ran_g[1],ran_g[-2]] ]
        #         if 'tixg' in kv_plt: tix_x = kv_plt['tixg']
        #         if len_e<=2 : tix_y = ran_e
        #         else        : tix_y = [y for y in [ran_e[1],ran_e[-2]] ]
        #         if 'tixe' in kv_plt: tix_y = kv_plt['tixe']
        #
        #         ''' % subplot y-axis % '''
        #         if midx==0:
        #             ax.tick_params(axis='y',direction='inout',left=True,labelleft=True,labelright=False,right=False)
        #             ax.yaxis.set_label_position('left')
        #             ax.set_ylabel(lab_y)
        #
        #             ax_pat = Patch(facecolor='none',edgecolor='none',label=labTex['leg'])
        #             ax.legend(handles=[ax_pat],loc='upper left',bbox_to_anchor=(kv_fig['patX'],kv_fig['patY']),frameon=True, fancybox=False,labelcolor=kv_fig['patC'],framealpha=0.0)
        #
        #         else: ax.tick_params(axis='y',direction='inout',left=True,labelleft=False,labelright=False,right=False)
        #         ax.set_yticks(tix_y), ax.set_yticklabels([f'{y-abs(evLo)}' for y in ax.get_yticks()])
        #         ''' % subplot x-axis % '''
        #
        #         if didx==self.dlen-1:
        #             ax.tick_params(axis='x',direction='in',top=False,labeltop=False,labelbottom=True,bottom=True)
        #             ax.xaxis.set_label_position('bottom'), ax.set_xticks(tix_x)
        #             ax.set_xlabel(lab_x), ax.set_xticklabels([f'{x-abs(gvLo)}' for x in ax.get_xticks() ])
        #
        #         elif didx!=0:
        #             ax.tick_params(axis='x',direction='in',top=False,labeltop=False,labelbottom=False,bottom=True)
        #             ax.xaxis.set_label_position('bottom'),ax.set_xticks(tix_x)
        #
        #         bar = ax.imshow(self.FCF_mtrx,norm=norm,origin='lower',cmap='plasma', aspect='auto')
        #         if self.kv_dd['didx']==self.dlen-1 and self.kv_dd['midx']==0:
        #             bar_h,bar_w = kv_fig['barH'],kv_fig['barW']
        #             bar_y,bar_x = kv_fig['barY'],kv_fig['barX']
        #             if midx==1: bar_x = kv_fig['barS']
        #             cbox = fig.add_axes([bar_x,bar_y,bar_w,bar_h])
        #             cbar = fig.colorbar(bar, cax=cbox, orientation=kv_fig['barO'])
        #             cbar.ax.text(kv_fig['barTx'],kv_fig['barTy'],labTex['cbar'],transform=cbar.ax.transAxes,ha='center',va='bottom',rotation=0,fontsize=12)
        #
        #     def fig_dim(self, *argv):
        #         kv_plt = self.kv_plt
        #         # Determine vibrational dimensions
        #         dim_ve = 58 if isinstance(kv_plt['ve'][0], str) else len(np.arange(min(kv_plt['ve']), max(kv_plt['ve']) + 1))
        #         dim_vg = 64 if isinstance(kv_plt['ψg'][0], str) else len(np.arange(min(kv_plt['ψg']), max(kv_plt['ψg']) + 1))
        #         rat_eg = dim_ve / dim_vg
        #         rat_ge = dim_vg / dim_ve
        #         # Grid shape
        #         rowNr = self.dlen ; colNr = self.mlen
        #         # Parametric scaling
        #         base_h_per_vg = 0.1
        #         base_w_per_ve = 0.1
        #         patC = clr_g
        #
        #         if kv_plt['fcf'] and kv_plt['fcf']=='plt':
        #             base_w_per_ve = 0.1
        #             base_h_per_vg = 0.25
        #             panel_w = 1   # one ground state column
        #             panel_h = 4   # one excited state row
        #             sepX = 0.02; sepY = 0.04
        #             # barX = 0.135; barW = 0.50; barY =-0.175;
        #             barX = .65; barW = .005; barTx = 3.000; barO = 'vertical'
        #             barY = .11; barH = .765; barTy = 1.005; barS = barX+barW+0.0245
        #             patX =-0.10; patY = 1.1
        #
        #             supX =0.385; supY = 0.0
        #             ttlX =0.5  ; ttlY = 1.0
        #
        #             figH = base_h_per_vg * dim_ve * rowNr + 4
        #             figW = base_w_per_ve * dim_vg * colNr + 1
        #             rtoW = [dim_vg] * colNr + [dim_vg]
        #             rtoH = [.2*figW+dim_vg]*rowNr
        #
        #         if kv_plt['fcf'] and kv_plt['fcf']=='img':
        #             if dim_ve==1:
        #                 panel_w = 1.8   # one ground state column
        #                 panel_h = 1.8   # one excited state row
        #                 sepX = .01; sepY = 0.01
        #                 barX = .65; barW = .005; barTx = 3.000; barO = 'vertical'
        #                 barY = .11; barH = .765; barTy = 1.005; barS = barX+barW+0.0245
        #                 patX =-.50; patY = 1.3
        #
        #                 supX =0.385; supY = 0.0
        #                 ttlX =0.5  ; ttlY = 1.0
        #
        #                 figH = base_h_per_vg * dim_ve * rowNr + 1.5
        #                 figW = base_w_per_ve * dim_vg * colNr + 4.15
        #                 rtoW = [dim_vg] * colNr + [dim_vg]
        #                 rtoH = [.2*figW+dim_vg]*rowNr
        #
        #             else:
        #                 panel_w = 3.0   # one ground state column
        #                 panel_h = 1.7   # one excited state row
        #                 sepX = 0.03 ; sepY = 0.00
        #                 barX = 0.82; barW = 0.02; barTx = 1.800; barO = 'vertical' # 'vertical'
        #                 barY = 0.11; barH = 0.77; barTy = 1.005; barS = barX+barW+0.3
        #                 patX =-0.85; patY = 0.65
        #                 supX = 0.47; supY = 0.09 #0.92
        #                 ttlX = 0.50; ttlY = 1.00 #-5.5
        #
        #                 figW = panel_w * colNr + 1.5  # space for colorbar
        #                 figH = panel_h * rowNr + 1.0  # space for labels
        #                 rtoW = [1] * colNr + [0.25]  # narrow colorbar
        #                 rtoH = [1] * rowNr
        #
        #         return {'barTx':barTx,'barTy':barTy,'barS':barS,'barO':barO,'ttlX':ttlX,'ttlY':ttlY,'supX':supX,'supY':supY,'patC':patC,'patX':patX,'patY':patY
        #                 ,'barX':barX,'barY':barY,'barH':barH,'barW':barW,'rtoW':rtoW,'rtoH':rtoH,'rowNr':rowNr,'colNr':colNr,'figH': figH,'figW':figW
        #                 ,'sepX': sepX,'sepY': sepY}
        #
        # kv_dd = {'s12':'s0'}
        # arr_ve = [[-2,-2]]
        # for ve in arr_ve:
        #     plot_dix = {'fcf':'plt','ψg':[-64,-2],'ve':ve}
        #     data_analys('png',kv_plt=plot_dix,kv_inp=kv_dd)
