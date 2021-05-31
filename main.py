# This program creates a small GUI tool to read in the COSMO INPUT namelist files and calculate for the user the
# the size of their data output for the model run.  Written 1/5/21 by Katie Osterried

import tkinter
from tkinter import ttk, StringVar
from tkinter.filedialog import askdirectory


# TODO: Error handling: there is a little bit of handling of user inputs, but could definitely be expanded
# TODO: Add restart files: right now it only calculates the output files, not the restart files
# TODO: Add ASCII files? (YU*, M_*, etc ..): these are not very big, so might not be too important to add
# TODO: Add assimilation files? (Feedobs): right now it doesn't include any of the data assimilation output
# TODO: Support output on subdomains?: not yet coded, not sure if necessary

class GUI:
    def __init__(self, window):
        self.input_dirtext = StringVar()
        self.input_dir = ''

        window.title("COSMO output size calculator")
        window.geometry("700x300")

        ttk.Button(window, text="Path to INPUT_* files", width=20, command=lambda: self.set_input_paths()). \
            grid(row=0, ipadx=5, ipady=15)
        ttk.Entry(window, textvariable=self.input_dirtext, width=60).grid(row=0, column=1, ipadx=1, ipady=1)

        label = ttk.Label(window)
        label.grid(row=2, column=1, pady=10)
        label.configure(background='white')
        ttk.Button(window, text="Calculate output size", command=lambda: self.buttonclick(label)).grid(row=1, column=1,
                                                                                                       ipadx=5,
                                                                                                       ipady=15)

    def set_input_paths(self):
        self.input_dir = askdirectory()
        self.input_dirtext.set(self.input_dir)

    def buttonclick(self, label):
        cmr = CosmoModelRun()
        orgpath = self.input_dir + "/INPUT_ORG"
        iopath = self.input_dir + "/INPUT_IO"
        if not orgpath:
            label.configure(foreground='red')
            label['text'] = "Please enter a valid path to a INPUT_ORG file"
        elif not iopath:
            label.configure(foreground='red')
            label['text'] = 'Please enter a valid path to a INPUT_IO file'
        else:
            label.configure(foreground='black')
            memtot, sfx = cmr.get_memory_size(self.input_dir)
            label['text'] = 'COSMO model output size: ' + str(round(memtot, 2)) + ' ' + sfx


class CosmoModelRun:

    @classmethod
    def get_param(cls, filename, param, ignore_comments=True, occurrence=1, startind=None, endind=None):
        import re

        namelist_pattern = re.compile(
            r""" ( (?P<varname> [a-zA-Z]\w*)[ ]* = [ ]*         # this reads the variable name part, it has to start
                                                                # with a letter and can have one single space before '='
              (?P<arg>                                          # defines the general argument list
               (?P<tal> (([ ]? , [ ]?)? (['].*?[']))+ )  |      # defines the text argument list
               (?P<nal> (([ ]? , [ ]?)? ([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?))+)  | # defines the numerical argument
               (?P<bla> [.] [a-zA-Z]+ [.]) )                    # defines the bool-like argument
              [ ]? ,?) """, re.VERBOSE)

        """retrieve a parameter from a Fortran namelist file"""

        comt = '!'  # the comment character
        zcount = 0
        try:
            data = open(filename).readlines()
            data = data[startind:endind]
            for ind, line in enumerate(data):
                # remove comment if required
                if ignore_comments:
                    i_cmt = line.find(comt)
                    if i_cmt > -1: line = line[0:i_cmt]
                else:
                    line = line.replace(comt, ' ')

                # search for  patterns, i.e. assignments
                matchobj = re.finditer(namelist_pattern, line)
                if matchobj is None:
                    continue
                for assignment in matchobj:
                    if assignment.group('varname') == param:
                        zcount += 1
                        if zcount == occurrence:
                            if (not param == 'yvarml') and (not param == 'yvarpl') and (not param == 'yvarzl'):
                                return assignment.group('arg')
                            else:
                                out_param = assignment.group('arg') + ','
                                inext = ind + 1
                                lnext = data[inext]
                                while '=' not in lnext:
                                    out_param += lnext
                                    inext = inext + 1
                                    lnext = data[inext]
                                return out_param

            return ''
        except FileNotFoundError:
            print("File not found:" + filename)

    def get_constant_vars(self, input_dir):
        nc2d = 7
        nc3d = 1
        # The following code is based on lines 5807 - 5895 of organize_data.f90 (git commit 4c27baa)
        lbdclim = self.get_param(input_dir + "/INPUT_IO", 'lbdclim')
        if not lbdclim or lbdclim.lower() == ".false.":
            nc2d = nc2d + 5
        llake = self.get_param(input_dir + "/INPUT_PHY", 'llake')
        if llake.lower() == ".true.":
            nc2d = nc2d + 2
        lforest = self.get_param(input_dir + "/INPUT_PHY", 'lforest')
        if lforest.lower() == ".true.":
            nc2d = nc2d + 2
        itype_canopy = self.get_param(input_dir + "/INPUT_PHY", 'itype_canopy')
        cskinc = self.get_param(input_dir + "/INPUT_PHY", 'cskinc')
        if int(itype_canopy == 2) and float(cskinc) < 0:
            nc2d = nc2d + 1
        lradtopo = self.get_param(input_dir + "/INPUT_PHY", 'lradtopo')
        if lradtopo.lower() == ".true.":
            nc2d = nc2d + 3
            nc3d = nc3d + 1
        itype_hydmod = self.get_param(input_dir + "/INPUT_PHY", 'itype_hydmod')
        if itype_hydmod and int(itype_hydmod) == 1:
            nc2d = nc2d + 1
        lsso = self.get_param(input_dir + "/INPUT_PHY", 'lsso')
        itype_vdif = self.get_param(input_dir + "/INPUT_PHY", 'itype_vdif')
        if lsso.lower() == ".true.":
            nc2d = nc2d + 4
        elif itype_vdif and int(itype_vdif) > -2:
            nc2d = nc2d + 1
        lemiss = self.get_param(input_dir + "/INPUT_PHY", 'lemiss')
        if lemiss.lower() == ".true.":
            nc2d = nc2d + 1
        lstomata = self.get_param(input_dir + "/INPUT_PHY", 'lstomata')
        if lstomata.lower() == ".true.":
            nc2d = nc2d + 1
        itype_aerosol = self.get_param(input_dir + "/INPUT_PHY", 'itype_aerosol')
        if itype_aerosol and int(itype_aerosol) == 2:
            nc2d = nc2d + 5
        elif itype_aerosol and int(itype_aerosol) == 3:
            nc3d = nc3d + 3
        itype_albedo = self.get_param(input_dir + "/INPUT_PHY", 'itype_albedo')
        if itype_albedo and int(itype_albedo) == 2:
            nc2d = nc2d + 2
        elif itype_albedo and int(itype_albedo) == 3:
            nc2d = nc2d + 1

        return nc2d, nc3d

    def get_gribout_indices(self, input_dir):
        iopath = input_dir + "/INPUT_IO"
        f_contents = open(iopath, 'r').read()
        numbers = []
        for num, line in enumerate(f_contents.strip().split('\n')):
            if "GRIBOUT" in line:
                numbers.append(num)
        return numbers

    def get_memory_size(self, input_dir):
        # Get the number of Gribout sections
        orgpath = input_dir + "/INPUT_ORG"
        iopath = input_dir + "/INPUT_IO"
        ngribout = self.get_param(iopath, 'ngribout')
        nx = int(self.get_param(orgpath, 'ie_tot'))
        ny = int(self.get_param(orgpath, 'je_tot'))
        nz = int(self.get_param(orgpath, 'ke_tot'))
        gmemtot = 0

        # Get the indices of the beginning of the gribout sections
        inds = self.get_gribout_indices(input_dir)
        # For each gribout section
        for x in range(int(ngribout)):
            gribout = GribOut()
            # Calculate the memory size
            if x < int(ngribout) - 1:
                gmem = gribout.get_gribout_memory_size(input_dir, inds[x], inds[x + 1], nx, ny, nz)
            else:
                gmem = gribout.get_gribout_memory_size(input_dir, inds[x], None, nx, ny, nz)

            # Sum up the total memory size
            gmemtot = gmemtot + gmem

        # Account for constant file
        nc2d, nc3d = self.get_constant_vars(input_dir)
        gmemtot = gmemtot + nc2d * nx * ny * 8 + nc3d * nx * ny * nz

        if gmemtot / 1000000 < 1:
            sfx = 'KB'
            gmemtot = gmemtot / 1000
        elif gmemtot / 1000000 > 1 > gmemtot / 1000000000:
            sfx = 'MB'
            gmemtot = gmemtot / 1000000
        elif gmemtot / 1000000000 > 1 > gmemtot / 1000000000000:
            sfx = 'GB'
            gmemtot = gmemtot / 1000000000
        elif gmemtot / 1000000000000 > 1 > gmemtot / 1000000000000000:
            sfx = 'TB'
            gmemtot = gmemtot / 1000000000000
        else:
            sfx = 'Bytes'
        return gmemtot, sfx


class GribOut:

    def get_gribout_memory_size(self, input_dir, startind, endind, nx, ny, nz):
        orgpath = input_dir + "/INPUT_ORG"
        iopath = input_dir + "/INPUT_IO"
        phypath = input_dir + "/INPUT_PHY"

        nout = self.get_output_time_steps(orgpath, iopath)
        n2dml, n3dml, lvarml = self.get_num_var(iopath, startind, endind, 'yvarml')

        n2dpl, n3dpl, lvarpl = self.get_num_var(iopath, startind, endind, 'yvarpl')
        pl = CosmoModelRun.get_param(iopath, 'plev', startind=startind, endind=endind)
        pl = len(pl.split(","))

        n2dzl, n3dzl, lvarzl = self.get_num_var(iopath, startind, endind, 'yvarzl')
        zl = CosmoModelRun.get_param(iopath, 'zlev', startind=startind, endind=endind)
        zl = len(zl.split(","))

        # TODO: Account for file headers
        # netcdf_io.f90: starting line: 3051
        # ie x je: zlatitude, zlongitude
        # ie: zrotlat
        # je: zrotlon
        # zcoord length: zvcoord
        # nwdirsec: zwdirsec_vals (nwdirsec hardwired to 16)
        # nbnds x nwdirsec: zwdirsec_bnds (nbnds hardwired to 2)
        # ntime: time (ntime hardwired to 1)
        # nbnds x ntime: time_bnds
        # If not luvmasspoint(?) and not constant (luvmasspoint in INPUT_IO)
        # ie x je: zslonu, zslatu, zslonv, zslatv
        # ie: zsrlon
        # je: zsrlat
        # If not constant file and not p or z file:
        # zheight_2m, zheight_10m, zheight_toa, zwbtemp_13c
        # ke_soil+1:czmls (ke_soil in INPUT_PHY)
        # nbnds x ke_soil+1: zsoil_bnds
        # If lmulti_snow
        # ke_snow: zind_snow (ke_snow in INPUT_PHY)
        # If RADARFWO and not p or z or s file and (ldo_bubbles or ldo_composite)
        # ktop-kbot+1: slev (ktop dependent on ldo_* flags, kbot = 1)
        # If lradtopo and not p or z
        # nhori: zsect (nhori in INPUT_PHY)
        # If lcalc_echotopz or lcalc_echotopp and not p and not z and not s
        # nechotop: dbzthresh_echotop (calculated in organize_data line 4662)
        # If luse_rttov and s
        # nmsgchan: zmsgchan_wave

        nvarml = 0
        if lvarml:
            nvarml = 1
        nvarpl = 0
        if lvarpl:
            nvarpl = 1
        nvarzl = 0
        if lvarzl:
            nvarzl = 1

        # First, account for global header variables
        glob_head_mem = ((nx * ny) * 8 + nx * 8 + ny * 8 + 16 * 8 + 16 * 2 * 8 + 8 + 16)(nvarml + nvarpl + nvarzl)
        luvmasspoint = CosmoModelRun.get_param(iopath, 'luvmasspoint')
        if luvmasspoint.lower() == '.true':
            glob_head_mem += (4 * (nx + ny) * 8 + nx * 8 + ny * 8)(nvarml + nvarpl + nvarzl)
            
        # Now, account for file specific headers
        head_mem = (nvarml*(nz + 4*8) + nvarpl*pl + nvarzl*zl) + glob_head_mem

        ke_soil = CosmoModelRun.get_param(phypath, 'ke_soil')
        if ke_soil:
            head_mem += (int(ke_soil) + 1 + 2*(int(ke_soil) + 1))*8

        lmulti_snow = CosmoModelRun.get_param(phypath, 'lmulti_snow')
        if lmulti_snow.lower() == ".true":
            ke_snow = CosmoModelRun.get_param(phypath, 'ke_snow')
            head_mem += int(ke_snow)

        mem_size = nout * (n2dml + n2dpl + n2dzl) * nx * ny * 8 + nout * (
                n3dml * nz + n3dpl * pl + n3dzl * zl) * nx * ny * 8
        return mem_size

    def get_num_var(self, iopath, startind, endind, var_type):
        import re
        num2d = 0
        num3d = 0
        lvar = False

        str2d = ['PMSL', 'DPSDT', 'FIS', 'HSURF', 'TO3', 'T_2M', 'T_2M_AV', 'TMAX_2M', 'TMIN_2M', 'TD_2M', 'TD_2M_AV',
                 'DD_10M', 'SP_10M', 'SP_10M_AV', 'U_10M', 'U_10M_AV', 'V_10M', 'V_10M_AV', 'QV_2M', 'RELHUM_2M', 'TQV',
                 'AEVAP_S', 'TQI', 'TOT_PR', 'TOT_PREC', 'CLCT', 'CLCT_AV', 'CLCL', 'CLCM', 'CLCH', 'TQC', 'SNOW_CON',
                 'SNOW_GSP', 'FR_LAND', 'Z0', 'ALB_RAD', 'T_CL', 'W_CL', 'PLCOV', 'RUNOFF_S', 'RUNOFF_G', 'FR_ICE',
                 'SNOW_MELT', 'SOBS_RAD', 'ASOB_S', 'THBS_RAD', 'ATHB_S', 'SOBT_RAD', 'ASOB_T', 'THBT_RAD', 'ATHB_T',
                 'LHFL_S', 'ALHFL_S', 'SHFL_S', 'ASHFL_S', 'UMFL_S', 'AUMFL_S', 'VMFL_S', 'AVMFL_S', 'PABS_RAD',
                 'APAB_S', 'ALHFL_BS', 'DURSUN', 'RSTOM', 'SWDIRS_RAD', 'ASWDIR_S', 'SWDIFDS_RAD', 'ASWDIFD_S',
                 'SWDIFUS_RAD', 'ASWDIFU_S', 'THDS_RAD', 'ATHD_S', 'THUS_RAD', 'ATHU_S', 'SODT_RAD', 'ASOD_T',
                 'SODS_RAD', 'ASOD_S', 'ASWDIR_SN', 'TQR', 'TQS', 'TQG', 'TWATER', 'TDIV_HUM', 'TCOND_MAX',
                 'TCOND10_MX', 'HBAS_SC', 'HTOP_SC', 'HBAS_CON', 'HTOP_CON', 'BAS_CON', 'TOP_CON', 'HTOP_DC', 'HZEROCL',
                 'SNOWLMT', 'GAMSO_LK', 'DP_BS_LK', 'DEPTH_LK', 'FETCH_LK', 'PRR_GSP', 'PRS_GSP', 'RAIN_GSP', 'PRR_CON',
                 'PRS_CON', 'RAIN_CON', 'TOT_SNOW', 'FRESHSNW', 'PRG_GSP', 'GRAU_GSP', 'PRH_GSP', 'HAIL_GSP', 'TQH',
                 'SDI_1', 'SDI_2', 'CAPE_MU', 'CIN_MU', 'CAPE_ML', 'CIN_ML', 'TKE_CON', 'TCM', 'TCH', 'VMAX_10M',
                 'T_BS_LK', 'LPI', 'LPI_MAX', 'SNOW_C', 'RSMIN', 'VABSMX_10M', 'VGUST_DYN_10M', 'VGUST_CON_10M',
                 'ASWDIR_SNO', 'DBZ_850', 'DBZ_CMAX', 'DBZ_CTMAX', 'EVATRA_SUM', 'TRA_SUM', 'TOTFORCE_S', 'RESID_WSO',
                 'MFLX_CON', 'CAPE_CON', 'MCONV', 'MCNV_CTMAX', 'QCVG_CON', 'FR_PAVED', 'AHF', 'SKC', 'SSO_STDH',
                 'SSO_GAMMA', 'SSO_THETA', 'SSO_SIGMA', 'FR_LAKE', 'EMIS_RAD', 'SOILTYP', 'LAI', 'ROOTDP', 'HMO3',
                 'VIO3', 'PLCOV_MX', 'PLCOV_MN', 'LAI_MX', 'LAI_MN', 'FOR_E', 'FOR_D', 'AER_SO4', 'AER_DUST', 'AER_ORG',
                 'AER_BC', 'AER_SS', 'TWATFLXU', 'ATWATFLXU', 'TWATFLXV', 'ATWATFLXV', 'SWDIR_COR', 'SLO_ANG',
                 'SLO_ASP', 'SKYVIEW', 'QVSFLX', 'FC', 'RLAT', 'RLON', 'ZTD', 'ZWD', 'ZHD', 'ALB_DRY', 'ALB_SAT',
                 'ALB_DIF', 'USTR_SSO', 'AUSTR_SSO', 'VSTR_SSO', 'AVSTR_SSO', 'VDIS_SSO', 'AVDIS_SSO', 'SOBS_RAD_CS',
                 'ASOB_S_CS', 'SOBT_RAD_CS', 'ASOB_T_CS', 'THBS_RAD_CS', 'ATHB_S_CS', 'THBT_RAD_CS', 'ATHB_T_CS',
                 'SUN_EL', 'SUN_AZI', 'UH_MAX', 'VORW_CTMAX', 'W_CTMAX', 'WT_DEPTH', 'W_SO_SL', 'LCL_ML', 'LFC_ML',
                 'CAPE_3KM', 'SWISS00', 'SWISS12', 'SLI', 'SI', 'BRN', 'HPBL', 'CEILING', 'CLDEPTH', 'CLCT_MOD',
                 'PMSL_ANAI', 'TQV_ANAI', 'TQC_ANAI', 'S_ORO_MAX', 'DURSUN_M', 'DURSUN_R', 'RAPA_SPPT', 'DHAIL_AV',
                 'DHAIL_AV', 'DHAIL_MX', 'W_UP_DUR', 'W_UP_MASK']

        yvar = CosmoModelRun.get_param(iopath, var_type, startind=startind, endind=endind)
        if yvar:
            yvar = yvar.split(",")
            lvar = True
            for outvar in yvar:
                var = re.sub('[^A-Z0-9_]', '', outvar)
                if var in str2d:
                    num2d = num2d + 1
                else:
                    num3d = num3d + 1

        return num2d, num3d, lvar

    def get_output_time_steps(self, orgpath, iopath):
        out_steps = []
        dt = CosmoModelRun.get_param(orgpath, 'dt')
        # Determine the last time step of the model run
        nstop = CosmoModelRun.get_param(orgpath, 'nstop')
        if not nstop:
            hstop = float(CosmoModelRun.get_param(orgpath, 'hstop'))
            nstop = int(round(float(hstop * 3600 / float(dt))))
        else:
            hstop = int(nstop) * float(dt) / 3600

        ngrib = CosmoModelRun.get_param(iopath, 'ngrib')
        if ngrib:
            out_steps.append(ngrib)

        hgrib = CosmoModelRun.get_param(iopath, 'hgrib')
        if hgrib:
            out_steps.append(int(round(float(hgrib) * 3600 / float(dt))))

        ncomb = CosmoModelRun.get_param(iopath, 'ncomb')
        if ncomb:
            ncomb = ncomb.split(',')
            nend = min(float(nstop), float(ncomb[1]))
            nnow = float(ncomb[0])
            while nnow < nend:
                out_steps.append(nnow)
                nnow = nnow + float(ncomb[2])

        hcomb = CosmoModelRun.get_param(iopath, 'hcomb')
        if hcomb:
            hcomb = hcomb.split(',')
            hend = min(float(hstop), float(hcomb[1]))
            hnow = float(hcomb[0])
            while hnow < hend:
                out_steps.append(int(round(float(hnow) * 3600 / float(dt))))
                hnow = hnow + float(hcomb[2])

        # Sort the out_steps
        GribOut.bubblesort(out_steps)

        # Eliminate duplicate entries
        list(dict.fromkeys(out_steps))
        # Calculate the number of output steps
        num_out = len(out_steps)
        return num_out

    @staticmethod
    def bubblesort(arr):
        n = len(arr)

        # Traverse through all array elements
        for i in range(n):

            # Last i elements are already in place
            for j in range(0, n - i - 1):

                # traverse the array from 0 to n-i-1
                # Swap if the element found is greater
                # than the next element
                if arr[j] > arr[j + 1]:
                    arr[j], arr[j + 1] = arr[j + 1], arr[j]


if __name__ == '__main__':
    window = tkinter.Tk()
    style = ttk.Style()
    style.theme_use('classic')
    gui = GUI(window)
    window.mainloop()
