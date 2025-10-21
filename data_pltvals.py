from dsty_str import TEX_gpsi,TEX_epsi,TEX_vB,TEX_c6,TEX_eΩ
from dsty_str import STY_eΩ,STY_vB,STY_ss,STY_Dy,STY_eΩ
from dsty_str import STYsí,STYno,STYrarr,STYlarr,STYdarr,STYsp,STYud,STYdt,STYds,STYerr,STYdos,STYsrc,STYdate,STYrun,STYmsg
# ***************************************************************************************************************************
from data_getvals import MFG_out
# ***************************************************************************************************************************
import os
path_loc = '.' ; fold_loc = os.listdir(path_loc)
# ***************************************************************************************************************************
import numpy as np
import math
import time
import seaborn as sns
c_sns = sns.color_palette("Paired",10)
import matplotlib.colors as colors
import matplotlib           ; import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.patches import Patch
from matplotlib.backends.backend_pdf import PdfPages
from IPython.display import display, Latex
from matplotlib.table import Table
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
# ***************************************************************************************************************************
plt.style.use('classic') , matplotlib.rcParams.update(matplotlib.rcParamsDefault) ; mpl.rcParams['savefig.dpi'] = 600
fsz = 11.5 ; lWd = .5 ; mSz = 3
plt.rcParams.update({'text.usetex':True ,'text.latex.preamble': r'\usepackage{newtxtext,newtxmath}\renewcommand{\familydefault}{\sfdefault}'
                    ,'figure.titlesize':fsz,'font.size':fsz-1
                    ,'legend.fontsize':fsz-4,'axes.titlesize':fsz-3
                    ,'ytick.labelsize':fsz-7,'xtick.labelsize':fsz-7
                    # ,'xtick.direction':'inout','ytick.direction':'out'
                    ,'savefig.transparent':True,'savefig.facecolor':'0.8'
                    ,'figure.max_open_warning':50,'figure.constrained_layout.use':True
                    })

clr_g = 'lightslategrey'


class MFG_plt():
    def __init__(self,*argv,run_str,**kwargs):

        iso_arr = [162,164]
        if not argv:

            for mm in iso_arr:
                dim = MFG_out('dim',**{'iso':mm,'run':run_str})

        if argv:

            self.src = 'data_pltvals.py | class MFG_plt '
            self.argv = argv  ; self.kv = kwargs.copy()
            self.kv.update({'iso':[162,164],'run':run_str})
            self.dlen, self.mlen = 6, len(self.kv['iso']) ;
            self.run = run_str ; kk_codev,kk_coder = 'v2','run%s'%(self.kv['run'])
            STYrun  = '%srun%s'%(kk_codev,self.kv['run'])
            self.TEX_run = r'$\rm{%s}%s\rm{%s}$'%(kk_codev,STY_vB,kk_coder)

            ls_ge = []

            def_src='%s %s | def __init__ '%(STYsrc,self.src)
            prt_err='%s %s '%(STYerr,STYrun)
            for mm in [162,164]:
                dim = MFG_out('dim',**{'iso':mm,'run':run_str}).kv_init
                if not dim:
                    print(def_src), print(prt_err+' D.N.E.')
                if dim:
                    # print(def_src)
                    ge_str = []
                    kv_ge = {'ψdim_g':dim['gdim'],'ψdim_e':dim['edim']}
                    for ge in ['e','g']:
                        s_vge = ge+'(-N,-1)' if 'e' in ge else ge+'(-M,-1)'
                        k_vge = 'ψ'+ge ; k_vxy = 'y'+ge if 'e' in ge else 'x'+ge
                        k_dim='%sdim_%s'%(k_vge[0],k_vge[-1]) ; k_idx='%sidx_%s'%(k_vge[0],k_vge[-1])
                        k_len='%slen_%s'%(k_vge[0],k_vge[-1]) ; k_ran='%sran_%s'%(k_vge[0],k_vge[-1])
                        ls_vran,ls_vidx,ls_vlen = [],[],[]
                        for idx,vdim in enumerate(kv_ge[k_dim]):
                            vran = range(0,vdim)
                            if k_vge in kwargs and isinstance(kwargs[k_vge],str):
                                if idx==0: ge_str.append('%s'%(s_vge))
                            if k_vge in kwargs and not isinstance(kwargs[k_vge],str):
                                if idx==0: ge_str.append('%s[%s,%s]'%( ge,min(kwargs[k_vge]),max(kwargs[k_vge]) ))
                                vv=kwargs[k_vge] ; vran=range(vdim-abs(min(vv)),vdim-abs(max(vv))+1)
                            vmin = min(vran)-vdim
                            if max(vran)-vdim==-1: vmax=None
                            else: vmax=max(vran)-vdim+1
                            vidx=(vmin,vmax) ; ls_vran.append(vran), ls_vlen.append(len(vran)), ls_vidx.append(vidx)
                        kv_ge.update({k_idx:ls_vidx,k_ran:ls_vran,k_len:ls_vlen})
                    ls_ge.append(kv_ge)
            if ls_ge:
                self.ls_ψm  = ls_ge.copy()
                self.kv_fig = self.dim_fig()
                self.kk_plt = '%s %s[plt-FCF]%s'%(STYdate,STYrun,'_'.join(ge_str))
                # self.get_fig()
                print(argv)


    def ls_data(self):
        ls_data = []
        for mm in self.kv['iso']: kv_data=MFG_out(**{'iso':mm,'run':self.run}).data_kv ; ls_data.append(kv_data)
        if ls_data: return ls_data

    def get_ax(self,**kwargs):
        mlen = self.mlen; midx = self.midx; fig = self.fig
        dlen = self.dlen; didx = self.didx; fgs = self.fgs
        ax = fig.add_subplot(fgs[midx,didx])
        if didx==dlen-1:
            if 'ylab' in kwargs: ax.set_ylabel(kwargs['ylab'])

        if midx==0: ax.tick_params(axis='x',bottom=False,labelbottom=False,labeltop=False,top=False)
        else:
            ax.tick_params(axis='x',bottom=True,labelbottom=True,labeltop=False,top=False)
            if 'xlab' in kwargs: ax.set_xlabel(kwargs['xlab'])

        ax.tick_params(axis='y',left=False,labelleft=False,labelright=True,right=True), ax.yaxis.set_label_position('right')
        if 'ytix' in kwargs: ax.set_yticks(kwargs['ytix'])
        if 'xtix' in kwargs: ax.set_xticks(kwargs['xtix'])
        if 'etix' in kwargs: ax.set_yticklabels(kwargs['etix'])
        if 'gtix' in kwargs: ax.set_xticklabels(kwargs['gtix'])
        if 'mlab' in kwargs and  didx==0: ax.set_title(kwargs['mlab'],x=-.45,y=0.7)
        if didx==self.dlen-1:
            im=ax.imshow(self.FCFmat,norm=self.FCFnorm,origin='lower',cmap='plasma',aspect='auto')
            chegh,cwdth=0.25,0.01
            cxpos=0.09

            ybase = 0.525
            if midx==0: cypos =  ybase
            if midx==1: cypos =  ybase-.4
            # cypos = cypos if midx==1 else cypos+.3
            im_cax =fig.add_axes([cxpos, cypos, cwdth,chegh])
            im_cbar=fig.colorbar(im, cax=im_cax, orientation='vertical')
            if 'clab' in kwargs:
                im_cbar.set_label(kwargs['clab'],fontsize=7)
                im_cbar.ax.yaxis.set_label_position('left'), im_cbar.ax.yaxis.tick_right()
        if 'plab' in kwargs and midx==0:
            ''' legend patch '''
            box_ψe = (-0.29,1.35) ; clr_ψe = 'black'
            pat_ψe = Patch(facecolor='none', edgecolor='none', label=kwargs['plab'])
            leg_ψe = ax.legend( handles=[pat_ψe], loc='upper left', bbox_to_anchor=box_ψe,
                                    frameon=True, fancybox=False, handletextpad=0.001, borderpad=0.001, labelspacing=0.001,
                                    labelcolor=clr_ψe, framealpha=0.0)# ; ax_fig.add_artist(leg_ψe)
        #
        # FCF_log = np.log(np.clip(FCF_mtrx, ψmin, ψmax))  # Correct: use zoomed-in matrix
        # ax.imshow(FCF_mtrx,norm=norm,origin='lower',cmap='plasma', aspect='auto')
        # pass
        return ax

    def get_fig(self):
        kv_fig=self.kv_fig
        fig=plt.figure(figsize=(kv_fig['figW'], kv_fig['figH']), dpi=100)
        fgs=gridspec.GridSpec(kv_fig['rowNr'],kv_fig['colNr']+1,figure=fig#,height_ratios=kv_fig['rtoH'],width_ratios=kv_fig['rtoW']
        )
        fig.subplots_adjust(hspace=kv_fig['sepY'],wspace=kv_fig['sepX'])
        self.fig,self.fgs=fig,fgs
        for midx,kv_dat in enumerate(self.ls_data()):
            self.midx = midx; kv_ψ = self.ls_ψm[midx]
            self.iso = self.kv['iso'][midx]


            TEX_mDy = r'$ ^{%s%g}%s$'%(STY_ss,self.iso,STY_Dy)
            TEX_run = self.TEX_run

            mlab = r'%s+%s'%(TEX_mDy,TEX_mDy)

            kv_ax = {'mlab':mlab,'xlab':TEX_gpsi,'ylab':TEX_epsi}
            # kv_ax = {'xlab':TEX_gpsi,'ylab':TEX_epsi}
            for didx,kv_dat in enumerate(kv_dat):
                self.didx = didx

                Ω_g=kv_dat['Ωg']; c6_g,s6_g=kv_dat['c6g'],kv_dat['s6g']; c12_g,s12_g=kv_dat['c12g'],kv_dat['s12g']
                Ω_e=kv_dat['Ωe']; c6_e,s6_e=kv_dat['c6e'],kv_dat['s6e']; c12_e,s12_e=kv_dat['c12e'],kv_dat['s12e']

                ψidx_e=kv_ψ['ψidx_e'][didx]; ψiLo_e,ψiHi_e=ψidx_e[0],ψidx_e[1]
                ψidx_g=kv_ψ['ψidx_g'][didx]; ψiLo_g,ψiHi_g=ψidx_g[0],ψidx_g[1]
                ψdim_e=kv_ψ['ψdim_e'][didx]; ψran_e=kv_ψ['ψran_e'][didx]
                ψdim_g=kv_ψ['ψdim_g'][didx]; ψran_g=kv_ψ['ψran_g'][didx]

                xtix = [ψran_g[ii] for ii in [0,-2]]; gtix = [f'{x-ψdim_g}' for x in xtix]
                ytix = [ψran_e[ii] for ii in [0,-2]]; etix = [f'{y-ψdim_e}' for y in ytix]

                FCFmat=kv_dat['FCFmat'][ψiLo_e:ψiHi_e,ψiLo_g:ψiHi_g]            ; self.FCFmat=FCFmat
                FCFval=np.concatenate([mat.ravel() for mat in FCFmat])
                FCFpos=FCFval[FCFval > 0]                                       ; ψmin,ψmax = max(1e-6,FCFpos.min()),FCFpos.max()
                norm=colors.LogNorm(vmin=ψmin, vmax=ψmax)                       ; self.FCFnorm=norm

                e6lab = '%s(%.0f)'%(TEX_c6,c6_e)
                eΩlab = '%s(%.0f)'%(TEX_eΩ,Ω_e)

                plab = r'%s%s$_{%s=%.0f}$'%(e6lab,TEX_vB,STY_eΩ,Ω_e)

                clab = mlab+'\n'+r'$\rm{log(FCF)}$'  if 'mlab' not in kv_ax else r'$\rm{log(FCF)}$'
                kv_ax.update({'xtix':xtix,'ytix':ytix,'gtix':gtix,'etix':etix,'plab':plab,'clab':clab})
                ax=self.get_ax(**kv_ax)
                ax.imshow(self.FCFmat,norm=self.FCFnorm,origin='lower',cmap='plasma', aspect='auto')

        for arg in self.argv:
            fig.suptitle(TEX_run,x=0.085 ,y=.99,color='red')
            if arg == 'plt': plt.show()
            if arg == 'png' or 'pdf':
                plt_name = self.kk_plt+'.'+arg ; plt_fold = '_plts'
                print(plt_name)
                os.makedirs(plt_fold, exist_ok=True)  # creates it if needed
                if arg=='png': plt.savefig(os.path.join(plt_fold, plt_name),bbox_inches='tight',facecolor='white',edgecolor='white')
                if arg=='pdf':
                    fig.patch.set_visible(False)
                    plt.savefig(os.path.join(plt_fold,plt_name),bbox_inches='tight',transparent=True)


    def dim_fig(self):
        kv_ψge,kv_fig = self.ls_ψm[0], {}
        len_ve,len_vg = min(kv_ψge['ψlen_e']),  min(kv_ψge['ψlen_g'])

        rowNr = self.mlen
        colNr = self.dlen
        kv_fig.update({'colNr':colNr,'rowNr':rowNr})

        # # Parametric scaling
        l_per_vg =0.01
        l_per_ve =0.01
        s_per_vg=len_vg*l_per_vg*0.1
        s_per_ve=len_ve*l_per_ve*0.1
        figH = 1.75#l_per_vg * len_ve * rowNr + 0
        figW = 9#l_per_ve * len_vg * colNr + 0
        kv_fig.update({'figW':figW,'figH':figH})

        rtoW = [s_per_ve] * colNr
        rtoH = [s_per_vg] * rowNr + [1]
        # kv_fig.update({'rtoW':rtoW,'rtoH':rtoH})

        sepX = 0.3 ; sepY = 0.04
        kv_fig.update({'sepX':sepX,'sepY':sepY})
        return kv_fig


run_set = ('A','B','C','D','E','F')
# run_set = ('A')
[MFG_plt(run_str=run,**{'ψe':'(-3,-1)','ψg':'(-2,-1)'}) for run in run_set]
