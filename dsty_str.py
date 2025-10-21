STYsí ='✓'
STYno ='✘'
STYrarr = '→'
STYlarr = '←'
STYdarr = '↓'
STYsp = ' '*100
STYud = '_'*100
STYdt = '.'*100
STYds = '-'*100
STYerr = '%s%s  err:'%(STYsp[0:3],STYrarr)
STYdos = '%s%s  dos:'%(STYsp[0:3],STYlarr)
STYsrc = '\n'+'%s src:'%(STYud[0:5])
STYset = '   %s  set:'%(STYsí)
STYrun = '   %s  run:'%(STYsí)
STYinp = '   %s  inp:'%(STYsí)
STYdne = '   %s  dne:'%(STYno)
STYmsg = '   %s     :'%(STYrarr)

import glob ; import datetime
stamp=(datetime.datetime.today())
str_date=str(stamp)[2:4]+str(stamp)[5:7]+str(stamp)[8:10] ; str_time=str(stamp)[11:13]+str(stamp)[14:16]+str(stamp)[17:19]
STYdate = ' '+str_date+':'+str_time

STY_pr = r'\prime'              ; STY_cdot = r'\cdot'
STY_s  = r'\scriptstyle'        ; STY_Lpar = r'\left(' ; STY_times = r'%s\times'%(STY_s)
STY_ss = r'\scriptscriptstyle'  ; STY_Rpar = r'\right)'
STY_vB = r'|'                   ; TEX_vB   = r'$%s$'%(STY_vB)
STY_vdB = r'\|'

STY_K9 = r'^{%s5}\rm{K}_{%s9}'%(STY_ss,STY_ss)      ; STY_0    = r'_{%s0}'%(STY_ss)
STY_I8 = r'^{%s5}\rm{I}_{%s8}'%(STY_ss,STY_ss)      ; STY_Dy   = r'\rm{Dy}'

TEX_I8  = r'$%s%s%s$'%(STY_Lpar,STY_I8,STY_Rpar)    ; TEX_Dy   = r'$\rm{Dy}$'
TEX_K9  = r'$%s%s%s$'%(STY_Lpar,STY_K9,STY_Rpar)    ;

TEX_FCF = r'$%s\rm{FCF}$'%(STY_s)
TEX_logFCF = r'FCF (log)'

STY_ve, STY_vg = r'v_{\rm e}', r'v_{\rm g}'
TEX_ve, TEX_vg = r'$%s$'%(STY_ve), r'$%s$'%(STY_vg)

STY_psi = r'\psi'
TEX_psi = r'$%s$'%STY_psi

TEX_pve = r'$%s(%s)$'%(STY_psi,STY_ve)
TEX_pvg = r'$%s(%s)$'%(STY_psi,STY_vg)

STY_e  = r'{\rm e}'
STY_g  = r'{\rm g}'
STY_Ω  = r'\Omega'
STY_gΩ = r'%s_{%s %s}'%(STY_Ω,STY_ss,STY_g)
STY_eΩ = r'%s_{%s %s}'%(STY_Ω,STY_ss,STY_e)

TEX_gpsi = r'$%s_{%s \rm g}$'%(STY_psi,STY_ss)
TEX_epsi = r'$%s_{%s \rm e}$'%(STY_psi,STY_ss)

TEX_vpr  = r'$v^{%s}$'%(STY_pr)

TEX_Ω  = r'$%s$'%(STY_Ω)
TEX_gΩ = r'$%s_{%s g}$'%(STY_Ω,STY_ss)
TEX_eΩ = r'$%s_{%s e}$'%(STY_Ω,STY_ss)

STY_Vlj = r'V^{%s{\rm LJ}}'%(STY_ss)
TEX_Vm  = r'$\rm V_{%s m}$'%(STY_ss)
TEX_Rm  = r'$\rm R_{%s m}$'%(STY_ss)

STY_C12  = r'C_{%s 12}'%(STY_ss)    ;
STY_c6   = r'C_{%s 6}'%(STY_ss)     ; TEX_c6  = r'$%s$'%(STY_c6)

TEX_Om  = r'$%s$'%(STY_Ω)
STY_pOm = r'%s_{%s e}'%(STY_Ω,STY_ss) ;
