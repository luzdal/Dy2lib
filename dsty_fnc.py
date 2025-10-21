import numpy as np
import os
path_loc = '.' ; fold_loc = os.listdir(path_loc)
import time
# ***************************************************************************************************************************

# ***************************************************************************************************************************
# print('holi')
def fnc_hash(arry): import hashlib; return hashlib.md5(arry).hexdigest() if isinstance(arry, bytes) else hashlib.md5(arry.tobytes()).hexdigest()

from scipy.constants import physical_constants, h, k, hbar, c, g
# from scipy.constants import u, convert_temperature, c, h, g, hbar, k, atm, bar, torr, mmHg, N_A

def get_const(str_):
    return physical_constants[str_][0]
m_p, m_e      = get_const('proton mass'), get_const('electron mass')

# Eh        = physical_constants['Hartree energy'][0]
# Eh_to_Hz  = Eh / h
# Eh_to_Mhz = Eh_to_Hz/1e6

kB , pi, sqrt = k, np.pi, np.sqrt
amu, cm2_m2   = 1.66e-27, 1e4

Eh = physical_constants['Hartree energy'][0]       # ~4.3597447222071e-18 J
a0 = physical_constants['Bohr radius'][0]          # ~5.29177210903e-11 m
# ***************************************************************************************************************************

def fnc_idxVm(**kwargs): return np.argmin(kwargs['Vpts'])
def fnc_idxV0(**kwargs): return min(range(len(['Vpts'])), key=lambda i: abs(kwargs['Vpts'][i]))

def fnc_slam(s): return s**(-1/6)
def fnc_lam(**kwargs): return (kwargs['c12']/kwargs['c6'])**(1/6)
def fnc_mass_mu(**kwargs): return 1/( 1/kwargs['m1'] + 1/kwargs['m2'] )
def fnc_RvdW(**kwargs): return a0*( abs(kwargs['c6'])*2*fnc_mass_mu(**kwargs)*m_p/m_e )**(1/4)/2
def fnc_EvdW(**kwargs): return h*hbar/(2*pi*2*fnc_mass_mu(**kwargs)*m_p) / fnc_RvdW(**kwargs)**2
def get_VLJ(**kwargs): return (kwargs['c12'] * kwargs['rpts']**(-12)) - (kwargs['c6'] * kwargs['rpts']**(-6))

def fnc_LJpot(**kwargs): return -16*(1-(kwargs['lam']/kwargs['rpts'])**6)/kwargs['rpts']**6
def fnc_ψIntpol(**kwargs):
    ''' interpolate excited-state wavefunctions to match ground grid '''
    from scipy.interpolate import interp1d as interpol
    ψe = np.zeros( ( len(kwargs['rpts_g']),kwargs['vpts'].shape[1] ) )
    for ψe_ii in range(kwargs['vpts'].shape[1]):
        get_interp=interpol(kwargs['rpts'],kwargs['vpts'][:,ψe_ii],kind='cubic',bounds_error=False,fill_value=0.0)
        ψe[:,ψe_ii]=get_interp(kwargs['rpts_g'])
    return ψe
def fnc_ψNormal(*argv,**kwargs):
    ''' normalize each vibrational wavefunction to unit probability '''
    ψN = np.zeros_like(kwargs['vpts'])
    for ii in range(kwargs['vpts'].shape[1]):
        norm = np.trapz(kwargs['vpts'][:, ii]**2, kwargs['rpts'])
        ψN[:, ii] = kwargs['vpts'][:, ii] / np.sqrt(norm) if norm > 0 else 0
    return ψN
def get_ψNorm(*argv,**kwargs):
    '''normalize(ψ_g or ψ_e_interpltd)'''
    dix_ge ={kk:vv for kk,vv in kwargs.items() if 'rpts' or 'vpts' in kk}
    if argv:
        ge = [argvs for argvs in argv if argvs=='g' or argvs=='e'][0]
        if ge =='g': return fnc_ψNormal(ge,**dix_ge)
        if ge =='e':
            ''' update ψ_e with interpolated values'''
            dix_ge.update({kk:fnc_ψIntpol(**dix_ge) for kk,vv in dix_ge.items() if 'vpts' in kk})
            return fnc_ψNormal(ge,**dix_ge)
def get_FCF(**kwargs):
    ''' computes Franck–Condon matrix from normalized wavefunctions '''
    FCFmat = np.zeros((kwargs['ψNor_e'].shape[1], kwargs['ψNor_g'].shape[1]))
    for mm in range(kwargs['ψNor_e'].shape[1]):
        for nn in range(kwargs['ψNor_g'].shape[1]):
            overlap = np.trapz( kwargs['ψNor_e'][:, mm] * kwargs['ψNor_g'][:, nn], kwargs['rpts_g'] )
            # if overlap!=0: print('FCF[%g,%g]:%g'%(mm,nn,np.abs(overlap)**2))
            FCFmat[mm,nn] = np.abs(overlap)**2
    return FCFmat
def open_file(file_path):
    data = []
    with open(file_path, 'r') as f:
        for line in f: row=np.fromstring(line.strip(), sep=' ') ; data.append(row)
    return data

def file_c6(*argv,**kwargs):
    for kk,vv in kwargs.items():
        file_name = () ; c6_txt = {'g':'vals_C6I8I8.txt', 'e':'vals_C6I8K9.txt'}
        if 'e' in kk: file_name = c6_txt['e'] ; file_dix = {'17':0,'16':0,'15':0}
        if 'g' in kk: file_name = c6_txt['g'] ; file_dix = {'16':0}
        if file_name:
            file_path = os.path.join(path_loc, file_name)
            file_vals = open_file(file_path)
            [print( 'Ω:%g , c6:%s'%(Ω_c6[0],Ω_c6[1:]) ) for Ω_c6 in file_vals if 'prt' in argv]
            # file_dix.update({str(int(Ω_c6[0])):Ω_c6[1:] for Ω_c6 in file_vals if str(int(Ω_c6[0])) in file_dix})
            file_dix.update({str(int(Ω_c6[0])): [float(val) for val in Ω_c6[1:]]
                                for Ω_c6 in file_vals
                                if str(int(Ω_c6[0])) in file_dix})


            if 'g' in kk: return file_dix[str(vv)][0]
            if 'e' in kk: return file_dix[str(vv)]

def val_mass(*argv,**kwargs):
    if 'iso' in kwargs and kwargs['iso']==162: return 161.9268056
    if 'iso' in kwargs and kwargs['iso']==164: return 163.9291819

def fnc_c12(*argv,**kwargs):
    c6,Rvdw,lam = kwargs['c6'], kwargs['Rvdw'],[vv for kk,vv in kwargs.items() if 'λ' in kk][0]
    return c6 * (lam*Rvdw)**6

def fnc_lam(*argv,**kwargs):
    lam = (kwargs['c12']/kwargs['c6'])**(1/6)/(kwargs['Rvdw']/a0)
    # print(lam)
    return lam

def Evdw_Rvdw(*argv,**kwargs):
    # λg = 9.115999998921005E-2 ; λe = 9.664006254444615E-2
    kv_val = {}
    if 'iso' in kwargs:
        mm_iso = val_mass(**kwargs)
        kv_val = {'m1':mm_iso,'m2':mm_iso}
    if 'amu' in kwargs:
        kv_val = {'m1':kwargs['amu'],'m2':kwargs['amu']}

    kv_val.update({'c6':vv for kk,vv in kwargs.items() if 'c6' in kk})
    ge = ''
    for kk,vv in kwargs.items():
        if 'g' in kk: ge = 'g'
        if 'e' in kk: ge = 'e'
    return {ge+'Evdw':fnc_EvdW(**kv_val)/Eh,ge+'Rvdw':fnc_RvdW(**kv_val)/a0}

def file_rename(*argv,**kwargs):
    old_name = kwargs.get('oldname')
    new_name = kwargs.get('newname')
    if not old_name or not new_name:
        print("Usage: file_rename(**{oldname='old.txt', newname='new.txt'})") ; return
    old_path = os.path.join(path_loc, old_name)
    new_path = os.path.join(path_loc, new_name)
    # Check if old file exists
    if not os.path.exists(old_path):
        print(f"Error: '{old_name}' does not exist in '{path_loc}'")
        return
    try:
        os.rename(old_path, new_path)
        print(f"Renamed: '{old_name}' → '{new_name}'")
    except Exception as e:
        print(f"Error: {e}")
def file_path(**kwargs):
    return os.path.join(kwargs['filepath'], kwargs['filename'])

    #
    # new_name = kwargs['newname']
    # old_name = kwargs['oldname']
    #
    # os.rename(os.path.join(fold_loc,old_name), os.path.join(fold_loc,new_name))

    # print(fold_loc)
    # print(kwargs)
    # if argv: print('file_rename(**{\'oldname\': ,\'newname\': })')
    # if 'oldname' and 'newname' in kwargs:
    #     os.rename(os.path.join(fold_loc, kwargs['oldname']), os.path.join(fold_loc, kwargs['newname'])) ; time.sleep(1)


    # print(kwargs)



    # for kk,vv in kwargs.items():
    #     if 'c6' in kk:
    #         print(kk,vv)
        # if 'amu' in kk: kv_vals.update({'m1':vv,'m2':vv})
    # print(kv_vals)
    #     if 'c6' in kk:
    #         c6 = vv
    #         Rvdw =
    #         print(c6)

            # kv_vals.update({'c6':vv})

        # print(kk[0],vv)
    # if len(kv_vals)==3:
    #     Rvdw = fnc_RvdW(**kv_vals) ;
    #     kv_vals.update(**{'Rvdw':fnc_RvdW(**kv_vals)})
    #     # print('')
    #     print(kv_vals)

    # dix_s = {'s0':1.00,'s1':0.95,'s2':1.05,'s3':1.02}
    # if kwargs:
    #     # print(kwargs)
    #     c12 = 361854239
    #     ge  = [kk[-1] for kk in kwargs.keys() for jj in ['c6'] if jj in kk][0]
    #
    #     dix = {'m1':kwargs['amu'],'m2':kwargs['amu']} ; dix.update(**{ij:vv for kk,vv in kwargs.items()
    #                                                                         for ij in ['Ω','c6'] if ij in kk})
    #
    #     Rvdw = fnc_RvdW(**dix) ; dix.update(**{'Rvdw':Rvdw})
    #     Evdw = fnc_EvdW(**dix) ; dix.update(**{'Evdw':Evdw})
    #
    #     if ge=='g': dix_s = {'s0':1.00}
    #     if ge=='e': dix_s = {'s0':1.00,'s1':0.95,'s2':1.05,'s3':1.02}
    #
    #     for argvs in argv:
    #         if 's' in argvs:
    #             kk = argvs ; vv = dix_s[kk]
    #             vv_c12 = int(c12*vv) ; vv_lam = fnc_lam(**dix,**{'c12':vv_c12})
    #             dix_kv = {kk:'%.2f'%vv,'c12_%s'%kk:vv_c12,'λ%s'%kk[-1]:vv_lam}
    #             dix.update(**{kk:'%.2f'%vv,'c12_%s'%kk:vv_c12,'λ%s'%kk[-1]:vv_lam})
    #             return dix
    # else:
#     #     for arg in argv: return dix_s[arg]
#
# def data_labels(*argv,**kwargs):
#     str_mm = '' ; Udash = '_' ; dash = '-'
#     if 'iso' in kwargs:
#         if round(kwargs['iso'])==162: str_mm = 'm1m1'
#         if round(kwargs['iso'])==164: str_mm = 'm2m2'
#         if 'mm' in argv: return str_mm
#
#         if 'c6' and 'Om' and 'lam' in kwargs:
#             val_c6 = kwargs['c6']; str_Om,str_λe = Udash+kwargs['Om'],Udash+kwargs['lam']
#             for kk,vv in file_c6('e').items():
#                 for idx_c6 in enumerate(vv):
#                     idx,c6 = idx_c6[0],idx_c6[1] ; cidx=dash+'c'+str(idx)
#                     if c6==val_c6: return '[%s]'%(str_mm+str_Om+cidx+str_λe)







#
#
