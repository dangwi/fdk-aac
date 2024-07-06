/* -----------------------------------------------------------------------------
Software License for The Fraunhofer FDK AAC Codec Library for Android

© Copyright  1995 - 2018 Fraunhofer-Gesellschaft zur Förderung der angewandten
Forschung e.V. All rights reserved.

 1.    INTRODUCTION
The Fraunhofer FDK AAC Codec Library for Android ("FDK AAC Codec") is software
that implements the MPEG Advanced Audio Coding ("AAC") encoding and decoding
scheme for digital audio. This FDK AAC Codec software is intended to be used on
a wide variety of Android devices.

AAC's HE-AAC and HE-AAC v2 versions are regarded as today's most efficient
general perceptual audio codecs. AAC-ELD is considered the best-performing
full-bandwidth communications codec by independent studies and is widely
deployed. AAC has been standardized by ISO and IEC as part of the MPEG
specifications.

Patent licenses for necessary patent claims for the FDK AAC Codec (including
those of Fraunhofer) may be obtained through Via Licensing
(www.vialicensing.com) or through the respective patent owners individually for
the purpose of encoding or decoding bit streams in products that are compliant
with the ISO/IEC MPEG audio standards. Please note that most manufacturers of
Android devices already license these patent claims through Via Licensing or
directly from the patent owners, and therefore FDK AAC Codec software may
already be covered under those patent licenses when it is used for those
licensed purposes only.

Commercially-licensed AAC software libraries, including floating-point versions
with enhanced sound quality, are also available from Fraunhofer. Users are
encouraged to check the Fraunhofer website for additional applications
information and documentation.

2.    COPYRIGHT LICENSE

Redistribution and use in source and binary forms, with or without modification,
are permitted without payment of copyright license fees provided that you
satisfy the following conditions:

You must retain the complete text of this software license in redistributions of
the FDK AAC Codec or your modifications thereto in source code form.

You must retain the complete text of this software license in the documentation
and/or other materials provided with redistributions of the FDK AAC Codec or
your modifications thereto in binary form. You must make available free of
charge copies of the complete source code of the FDK AAC Codec and your
modifications thereto to recipients of copies in binary form.

The name of Fraunhofer may not be used to endorse or promote products derived
from this library without prior written permission.

You may not charge copyright license fees for anyone to use, copy or distribute
the FDK AAC Codec software or your modifications thereto.

Your modified versions of the FDK AAC Codec must carry prominent notices stating
that you changed the software and the date of any change. For modified versions
of the FDK AAC Codec, the term "Fraunhofer FDK AAC Codec Library for Android"
must be replaced by the term "Third-Party Modified Version of the Fraunhofer FDK
AAC Codec Library for Android."

3.    NO PATENT LICENSE

NO EXPRESS OR IMPLIED LICENSES TO ANY PATENT CLAIMS, including without
limitation the patents of Fraunhofer, ARE GRANTED BY THIS SOFTWARE LICENSE.
Fraunhofer provides no warranty of patent non-infringement with respect to this
software.

You may use this FDK AAC Codec software or modifications thereto only for
purposes that are authorized by appropriate patent licenses.

4.    DISCLAIMER

This FDK AAC Codec software is provided by Fraunhofer on behalf of the copyright
holders and contributors "AS IS" and WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES,
including but not limited to the implied warranties of merchantability and
fitness for a particular purpose. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE for any direct, indirect, incidental, special, exemplary,
or consequential damages, including but not limited to procurement of substitute
goods or services; loss of use, data, or profits, or business interruption,
however caused and on any theory of liability, whether in contract, strict
liability, or tort (including negligence), arising in any way out of the use of
this software, even if advised of the possibility of such damage.

5.    CONTACT INFORMATION

Fraunhofer Institute for Integrated Circuits IIS
Attention: Audio and Multimedia Departments - FDK AAC LL
Am Wolfsmantel 33
91058 Erlangen, Germany

www.iis.fraunhofer.de/amm
amm-info@iis.fraunhofer.de
----------------------------------------------------------------------------- */

/******************* Library for basic calculation routines ********************

   Author(s): Dangwi Xu @ sophgo

   Description: dct RISC-V assembler replacements.

*******************************************************************************/

#if defined(__riscv_xthead) && __riscv_v == 7000
#include <riscv_vector.h>

#ifndef FUNCTION_dct_IV
#define FUNCTION_dct_IV

void dct_IV(FIXP_DBL *pDat, int L, int *pDat_e) {
  int sin_step = 0;
  int M = L >> 1;

  const FIXP_WTP *twiddle;
  const FIXP_STP *sin_twiddle;

  FDK_ASSERT(L >= 4);

  FDK_ASSERT(L >= 4);

  dct_getTables(&twiddle, &sin_twiddle, &sin_step, L);

  {
    FIXP_DBL *RESTRICT pDat_0 = &pDat[0];
    FIXP_DBL *RESTRICT pDat_1 = &pDat[L - 2];
    int i, vl, vlmax;
    vuint32m1_t vDatidx, vDatlen;
    vint32m1_t  vaacu1, vaacu2, vaacu3, vaacu4;
    vint32m1_t  vmul1, vmul2, vresult;
    vint32m1_t  v32twiddlep0, v32twiddlep1;
    vint32m1_t  v32rep0, v32rep1, v32imp0, v32imp1;

    vlmax   = vsetvlmax_e32m1();
    vDatlen = vmv_v_x_u32m1((L - 2) << 2, vlmax);
    vDatidx = vid_v_u32m1(vlmax);
    vDatidx = vsll_vx_u32m1(vDatidx, 3, vlmax);
    vDatidx = vsub_vv_u32m1(vDatlen, vDatidx, vlmax);

    vwrite_csr(RVV_VXRM, 0b10);
    for (i = 0; i < M - 1; i += 2 * vl) {
      vl = vsetvl_e32m1((M - i) >> 1);

      vaacu1 = vlxw_v_i32m1(&pDat[1], vDatidx, vl);
      vaacu2 = vlse32_v_i32m1(&pDat[i + 0], 8, vl);
      vaacu3 = vlse32_v_i32m1(&pDat[i + 1], 8, vl);
      vaacu4 = vlxw_v_i32m1(&pDat[0], vDatidx, vl);

      v32twiddlep0 = vlse32_v_i32m1((int32_t*)&twiddle[i + 0], 8, vl);
      v32twiddlep1 = vlse32_v_i32m1((int32_t*)&twiddle[i + 1], 8, vl);

      v32rep0 = vsll_vx_i32m1(v32twiddlep0, 16, vl);
      v32rep1 = vsll_vx_i32m1(v32twiddlep1, 16, vl);
      v32imp0 = vand_vx_i32m1(v32twiddlep0, 0xffff0000, vl);
      v32imp1 = vand_vx_i32m1(v32twiddlep1, 0xffff0000, vl);

      vmul1 = vmulh_vv_i32m1(vaacu1, v32rep0, vl);
      vmul2 = vmulh_vv_i32m1(vaacu2, v32imp0, vl);
      vresult = vasub_vv_i32m1(vmul1, vmul2, vl);
      vsse32_v_i32m1(&pDat[i + 1], 8, vresult, vl);

      vmul1 = vmulh_vv_i32m1(vaacu2, v32rep0, vl);
      vmul2 = vmulh_vv_i32m1(vaacu1, v32imp0, vl);
      vresult = vaadd_vv_i32m1(vmul1, vmul2, vl);
      vsse32_v_i32m1(&pDat[i + 0], 8, vresult, vl);

      vmul1 = vmulh_vv_i32m1(vaacu4, v32rep1, vl);
      vmul2 = vmulh_vv_i32m1(vaacu3, v32imp1, vl);
      vresult = vasub_vv_i32m1(vmul1, vmul2, vl);
      vresult = vneg_v_i32m1(vresult, vl);
      vsxw_v_i32m1(&pDat[1], vDatidx, vresult, vl);

      vmul1 = vmulh_vv_i32m1(vaacu3, v32rep1, vl);
      vmul2 = vmulh_vv_i32m1(vaacu4, v32imp1, vl);
      vresult = vaadd_vv_i32m1(vmul1, vmul2, vl);
      vsxw_v_i32m1(&pDat[0], vDatidx, vresult, vl);

      vDatidx = vsub_vx_u32m1(vDatidx, vl * 8, vl);
    }

    if (M & 1) {
      FIXP_DBL accu1, accu2;
      int idx0 = M & 0xfffffffe;
      int idx1 = L - 2 - (M & 0xfffffffe);

      accu1 = pDat[idx1 + 1];
      accu2 = pDat[idx0 + 0];

      cplxMultDiv2(&accu1, &accu2, accu1, accu2, twiddle[M - 1]);

      pDat[idx0 + 0] = accu2 >> 1;
      pDat[idx0 + 1] = accu1 >> 1;
    }
  }

  fft(M, pDat, pDat_e);

  {
    FIXP_DBL *RESTRICT pDat_0 = &pDat[0];
    FIXP_DBL *RESTRICT pDat_1 = &pDat[L - 2];
    FIXP_DBL accu1, accu2, accu3, accu4;
    int idx, i;

    /* Sin and Cos values are 0.0f and 1.0f */
    accu1 = pDat_1[0];
    accu2 = pDat_1[1];

    pDat_1[1] = -pDat_0[1];

    /* 28 cycles for ARM926 */
    for (idx = sin_step, i = 1; i<(M + 1)>> 1; i++, idx += sin_step) {
      FIXP_STP twd = sin_twiddle[idx];
      cplxMult(&accu3, &accu4, accu1, accu2, twd);
      pDat_0[1] = accu3;
      pDat_1[0] = accu4;

      pDat_0 += 2;
      pDat_1 -= 2;

      cplxMult(&accu3, &accu4, pDat_0[1], pDat_0[0], twd);

      accu1 = pDat_1[0];
      accu2 = pDat_1[1];

      pDat_1[1] = -accu3;
      pDat_0[0] = accu4;
    }

    if ((M & 1) == 0) {
      /* Last Sin and Cos value pair are the same */
      accu1 = fMult(accu1, WTC(0x5a82799a));
      accu2 = fMult(accu2, WTC(0x5a82799a));

      pDat_1[0] = accu1 + accu2;
      pDat_0[1] = accu1 - accu2;
    }
  }

  /* Add twiddeling scale. */
  *pDat_e += 2;
}
#endif // FUNCTION_dct_IV

#endif // defined(__riscv_xthead) && __riscv_v == 7000