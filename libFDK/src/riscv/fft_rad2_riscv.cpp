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

   Description: dit_fft RISC-V assembler replacements.

*******************************************************************************/

#ifndef __FFT_RAD2_CPP__
#error \
    "Do not compile this file separately. It is included on demand from fft_rad2.cpp"
#endif

#if defined(__riscv_xthead) && __riscv_v == 7000
#include <riscv_vector.h>

#ifndef FUNCTION_dit_fft
#define FUNCTION_dit_fft

void dit_fft(FIXP_DBL *x, const INT ldn, const FIXP_STP *trigdata,
             const INT trigDataSize) {
  const INT n = 1 << ldn;
  INT trigstep, i, ldm;
  C_ALLOC_ALIGNED_CHECK(x);

  scramble(x, n);
  /*
   * 1+2 stage radix 4
   */
  {
    INT vl;
    vint32m1_t vx0, vx1, vx2, vx3, vx4, vx5, vx6, vx7;
    vint32m1_t va0, va1, va2, va3;

    vwrite_csr(RVV_VXRM, 0b10);
    for (i = 0; i < n * 2; i += vl * 8) {
      vl  = vsetvl_e32m1((n * 2 - i) >> 3);
      vx0 = vlse32_v_i32m1(&x[i + 0], 4 * 8, vl);
      vx1 = vlse32_v_i32m1(&x[i + 1], 4 * 8, vl);
      vx2 = vlse32_v_i32m1(&x[i + 2], 4 * 8, vl);
      vx3 = vlse32_v_i32m1(&x[i + 3], 4 * 8, vl);
      vx4 = vlse32_v_i32m1(&x[i + 4], 4 * 8, vl);
      vx5 = vlse32_v_i32m1(&x[i + 5], 4 * 8, vl);
      vx6 = vlse32_v_i32m1(&x[i + 6], 4 * 8, vl);
      vx7 = vlse32_v_i32m1(&x[i + 7], 4 * 8, vl);

      va0 = vaadd_vv_i32m1(vx0, vx2, vl);
      va1 = vaadd_vv_i32m1(vx4, vx6, vl);
      va2 = vaadd_vv_i32m1(vx1, vx3, vl);
      va3 = vaadd_vv_i32m1(vx5, vx7, vl);

      vx0 = vadd_vv_i32m1(va0, va1, vl);
      vx4 = vsub_vv_i32m1(va0, va1, vl);
      vx1 = vadd_vv_i32m1(va2, va3, vl);
      vx5 = vsub_vv_i32m1(va2, va3, vl);

      va0 = vsub_vv_i32m1(va0, vx2, vl);
      va1 = vsub_vv_i32m1(va1, vx6, vl);
      va2 = vsub_vv_i32m1(va2, vx3, vl);
      va3 = vsub_vv_i32m1(va3, vx7, vl);

      vx2 = vadd_vv_i32m1(va0, va3, vl);
      vx6 = vsub_vv_i32m1(va0, va3, vl);
      vx3 = vsub_vv_i32m1(va2, va1, vl);
      vx7 = vadd_vv_i32m1(va2, va1, vl);

      vsse32_v_i32m1(&x[i + 0], 4 * 8, vx0, vl);
      vsse32_v_i32m1(&x[i + 1], 4 * 8, vx1, vl);
      vsse32_v_i32m1(&x[i + 2], 4 * 8, vx2, vl);
      vsse32_v_i32m1(&x[i + 3], 4 * 8, vx3, vl);
      vsse32_v_i32m1(&x[i + 4], 4 * 8, vx4, vl);
      vsse32_v_i32m1(&x[i + 5], 4 * 8, vx5, vl);
      vsse32_v_i32m1(&x[i + 6], 4 * 8, vx6, vl);
      vsse32_v_i32m1(&x[i + 7], 4 * 8, vx7, vl);
    }
  }

  for (ldm = 3; ldm <= ldn; ++ldm) {
    INT m = (1 << ldm);
    INT mh = (m >> 1);
    INT j, r;

    trigstep = ((trigDataSize << 2) >> ldm);

    FDK_ASSERT(trigstep > 0);

    /* Do first iteration with c=1.0 and s=0.0 separately to avoid loosing to
       much precision. Beware: The impact on the overal FFT precision is rather
       large. */
    { /* block 1 */

      j = 0;

      for (r = 0; r < n; r += m) {
        INT t1 = (r + j) << 1;
        INT t2 = t1 + (mh << 1);
        FIXP_DBL vr, vi, ur, ui;

        // cplxMultDiv2(&vi, &vr, x[t2+1], x[t2], (FIXP_SGL)1.0, (FIXP_SGL)0.0);
        vi = x[t2 + 1] >> 1;
        vr = x[t2] >> 1;

        ur = x[t1] >> 1;
        ui = x[t1 + 1] >> 1;

        x[t1] = ur + vr;
        x[t1 + 1] = ui + vi;

        x[t2] = ur - vr;
        x[t2 + 1] = ui - vi;

        t1 += mh;
        t2 = t1 + (mh << 1);

        // cplxMultDiv2(&vr, &vi, x[t2+1], x[t2], (FIXP_SGL)1.0, (FIXP_SGL)0.0);
        vr = x[t2 + 1] >> 1;
        vi = x[t2] >> 1;

        ur = x[t1] >> 1;
        ui = x[t1 + 1] >> 1;

        x[t1] = ur + vr;
        x[t1 + 1] = ui - vi;

        x[t2] = ur - vr;
        x[t2 + 1] = ui + vi;
      }

    } /* end of  block 1 */

    for (j = 1; j < mh / 4; ++j) {
      FIXP_STP cs;

      cs = trigdata[j * trigstep];

      if (ldn - ldm <= 1) { /* Small FFTs */
        for (r = 0; r < n; r += m) {
          INT t1 = (r + j) << 1;
          INT t2 = t1 + (mh << 1);
          FIXP_DBL vr, vi, ur, ui;

          cplxMultDiv2(&vi, &vr, x[t2 + 1], x[t2], cs);

          ur = x[t1] >> 1;
          ui = x[t1 + 1] >> 1;

          x[t1] = ur + vr;
          x[t1 + 1] = ui + vi;

          x[t2] = ur - vr;
          x[t2 + 1] = ui - vi;

          t1 += mh;
          t2 = t1 + (mh << 1);

          cplxMultDiv2(&vr, &vi, x[t2 + 1], x[t2], cs);

          ur = x[t1] >> 1;
          ui = x[t1 + 1] >> 1;

          x[t1] = ur + vr;
          x[t1 + 1] = ui - vi;

          x[t2] = ur - vr;
          x[t2 + 1] = ui + vi;

          /* Same as above but for t1,t2 with j>mh/4 and thus cs swapped */
          t1 = (r + mh / 2 - j) << 1;
          t2 = t1 + (mh << 1);

          cplxMultDiv2(&vi, &vr, x[t2], x[t2 + 1], cs);

          ur = x[t1] >> 1;
          ui = x[t1 + 1] >> 1;

          x[t1] = ur + vr;
          x[t1 + 1] = ui - vi;

          x[t2] = ur - vr;
          x[t2 + 1] = ui + vi;

          t1 += mh;
          t2 = t1 + (mh << 1);

          cplxMultDiv2(&vr, &vi, x[t2], x[t2 + 1], cs);

          ur = x[t1] >> 1;
          ui = x[t1 + 1] >> 1;

          x[t1] = ur - vr;
          x[t1 + 1] = ui - vi;

          x[t2] = ur + vr;
          x[t2 + 1] = ui + vi;
        }
      }
      else { /* vector FFTs */
        int vl;
        FIXP_DBL csvre, csvim;
        vint32m1_t vxt1p0, vxt1p1;
        vint32m1_t vxt2p0, vxt2p1;
        vint32m1_t vmul1, vmul2;
        vint32m1_t vvi, vvr, vur, vui;

        csvre = FX_SGL2FX_DBL(cs.v.re);
        csvim = FX_SGL2FX_DBL(cs.v.im);
        for (r = 0; r < n; r += m * vl) {
          vl = vsetvl_e32m1((n - r) >> ldm);

          // FFT round1
          vxt1p0 = vlse32_v_i32m1(&x[((r + j) << 1) + 0], 8 * m, vl);
          vxt1p1 = vlse32_v_i32m1(&x[((r + j) << 1) + 1], 8 * m, vl);
          vxt2p0 = vlse32_v_i32m1(&x[((r + j + mh) << 1) + 0], 8 * m, vl);
          vxt2p1 = vlse32_v_i32m1(&x[((r + j + mh) << 1) + 1], 8 * m, vl);

          vmul1 = vmulh_vx_i32m1(vxt2p1, csvre, vl);
          vmul2 = vmulh_vx_i32m1(vxt2p0, csvim, vl);
          vvi = vsub_vv_i32m1(vmul1, vmul2, vl);

          vmul1 = vmulh_vx_i32m1(vxt2p0, csvre, vl);
          vmul2 = vmulh_vx_i32m1(vxt2p1, csvim, vl);
          vvr = vadd_vv_i32m1(vmul1, vmul2, vl);

          vur = vsra_vx_i32m1(vxt1p0, 1, vl);
          vui = vsra_vx_i32m1(vxt1p1, 1, vl);

          vxt1p0 = vadd_vv_i32m1(vur, vvr, vl);
          vxt1p1 = vadd_vv_i32m1(vui, vvi, vl);

          vxt2p0 = vsub_vv_i32m1(vur, vvr, vl);
          vxt2p1 = vsub_vv_i32m1(vui, vvi, vl);

          vsse32_v_i32m1(&x[((r + j) << 1) + 0], 8 * m, vxt1p0, vl);
          vsse32_v_i32m1(&x[((r + j) << 1) + 1], 8 * m, vxt1p1, vl);
          vsse32_v_i32m1(&x[((r + j + mh) << 1) + 0], 8 * m, vxt2p0, vl);
          vsse32_v_i32m1(&x[((r + j + mh) << 1) + 1], 8 * m, vxt2p1, vl);

          // FFT round2
          vxt1p0 = vlse32_v_i32m1(&x[((r + j) << 1) + 0 + mh], 8 * m, vl);
          vxt1p1 = vlse32_v_i32m1(&x[((r + j) << 1) + 1 + mh], 8 * m, vl);
          vxt2p0 = vlse32_v_i32m1(&x[((r + j + mh) << 1) + 0 + mh], 8 * m, vl);
          vxt2p1 = vlse32_v_i32m1(&x[((r + j + mh) << 1) + 1 + mh], 8 * m, vl);

          vmul1 = vmulh_vx_i32m1(vxt2p1, csvre, vl);
          vmul2 = vmulh_vx_i32m1(vxt2p0, csvim, vl);
          vvr = vsub_vv_i32m1(vmul1, vmul2, vl);

          vmul1 = vmulh_vx_i32m1(vxt2p0, csvre, vl);
          vmul2 = vmulh_vx_i32m1(vxt2p1, csvim, vl);
          vvi = vadd_vv_i32m1(vmul1, vmul2, vl);

          vur = vsra_vx_i32m1(vxt1p0, 1, vl);
          vui = vsra_vx_i32m1(vxt1p1, 1, vl);

          vxt1p0 = vadd_vv_i32m1(vur, vvr, vl);
          vxt1p1 = vsub_vv_i32m1(vui, vvi, vl);

          vxt2p0 = vsub_vv_i32m1(vur, vvr, vl);
          vxt2p1 = vadd_vv_i32m1(vui, vvi, vl);

          vsse32_v_i32m1(&x[((r + j) << 1) + 0 + mh], 8 * m, vxt1p0, vl);
          vsse32_v_i32m1(&x[((r + j) << 1) + 1 + mh], 8 * m, vxt1p1, vl);
          vsse32_v_i32m1(&x[((r + j + mh) << 1) + 0 + mh], 8 * m, vxt2p0, vl);
          vsse32_v_i32m1(&x[((r + j + mh) << 1) + 1 + mh], 8 * m, vxt2p1, vl);

          // FFT round3
          vxt1p0 = vlse32_v_i32m1(&x[((r + mh / 2 - j) << 1) + 0], 8 * m, vl);
          vxt1p1 = vlse32_v_i32m1(&x[((r + mh / 2 - j) << 1) + 1], 8 * m, vl);
          vxt2p0 = vlse32_v_i32m1(&x[((r + mh / 2 - j + mh) << 1) + 0], 8 * m, vl);
          vxt2p1 = vlse32_v_i32m1(&x[((r + mh / 2 - j + mh) << 1) + 1], 8 * m, vl);

          vmul1 = vmulh_vx_i32m1(vxt2p0, csvre, vl);
          vmul2 = vmulh_vx_i32m1(vxt2p1, csvim, vl);
          vvi = vsub_vv_i32m1(vmul1, vmul2, vl);

          vmul1 = vmulh_vx_i32m1(vxt2p1, csvre, vl);
          vmul2 = vmulh_vx_i32m1(vxt2p0, csvim, vl);
          vvr = vadd_vv_i32m1(vmul1, vmul2, vl);

          vur = vsra_vx_i32m1(vxt1p0, 1, vl);
          vui = vsra_vx_i32m1(vxt1p1, 1, vl);

          vxt1p0 = vadd_vv_i32m1(vur, vvr, vl);
          vxt1p1 = vsub_vv_i32m1(vui, vvi, vl);

          vxt2p0 = vsub_vv_i32m1(vur, vvr, vl);
          vxt2p1 = vadd_vv_i32m1(vui, vvi, vl);

          vsse32_v_i32m1(&x[((r + mh / 2 - j) << 1) + 0], 8 * m, vxt1p0, vl);
          vsse32_v_i32m1(&x[((r + mh / 2 - j) << 1) + 1], 8 * m, vxt1p1, vl);
          vsse32_v_i32m1(&x[((r + mh / 2 - j + mh) << 1) + 0], 8 * m, vxt2p0, vl);
          vsse32_v_i32m1(&x[((r + mh / 2 - j + mh) << 1) + 1], 8 * m, vxt2p1, vl);

          // FFT round4
          vxt1p0 = vlse32_v_i32m1(&x[((r + mh / 2 - j) << 1) + 0 + mh], 8 * m, vl);
          vxt1p1 = vlse32_v_i32m1(&x[((r + mh / 2 - j) << 1) + 1 + mh], 8 * m, vl);
          vxt2p0 = vlse32_v_i32m1(&x[((r + mh / 2 - j + mh) << 1) + 0 + mh], 8 * m, vl);
          vxt2p1 = vlse32_v_i32m1(&x[((r + mh / 2 - j + mh) << 1) + 1 + mh], 8 * m, vl);

          vmul1 = vmulh_vx_i32m1(vxt2p0, csvre, vl);
          vmul2 = vmulh_vx_i32m1(vxt2p1, csvim, vl);
          vvr = vsub_vv_i32m1(vmul1, vmul2, vl);

          vmul1 = vmulh_vx_i32m1(vxt2p1, csvre, vl);
          vmul2 = vmulh_vx_i32m1(vxt2p0, csvim, vl);
          vvi = vadd_vv_i32m1(vmul1, vmul2, vl);

          vur = vsra_vx_i32m1(vxt1p0, 1, vl);
          vui = vsra_vx_i32m1(vxt1p1, 1, vl);

          vxt1p0 = vsub_vv_i32m1(vur, vvr, vl);
          vxt1p1 = vsub_vv_i32m1(vui, vvi, vl);

          vxt2p0 = vadd_vv_i32m1(vur, vvr, vl);
          vxt2p1 = vadd_vv_i32m1(vui, vvi, vl);

          vsse32_v_i32m1(&x[((r + mh / 2 - j) << 1) + 0 + mh], 8 * m, vxt1p0, vl);
          vsse32_v_i32m1(&x[((r + mh / 2 - j) << 1) + 1 + mh], 8 * m, vxt1p1, vl);
          vsse32_v_i32m1(&x[((r + mh / 2 - j + mh) << 1) + 0 + mh], 8 * m, vxt2p0, vl);
          vsse32_v_i32m1(&x[((r + mh / 2 - j + mh) << 1) + 1 + mh], 8 * m, vxt2p1, vl);
        }
      }
    }

    { /* block 2 */
      j = mh / 4;

      if (ldn - ldm <= 1) { /* Small FFTs */
        for (r = 0; r < n; r += m) {
          INT t1 = (r + j) << 1;
          INT t2 = t1 + (mh << 1);
          FIXP_DBL vr, vi, ur, ui;

          cplxMultDiv2(&vi, &vr, x[t2 + 1], x[t2], STC(0x5a82799a),
                      STC(0x5a82799a));

          ur = x[t1] >> 1;
          ui = x[t1 + 1] >> 1;

          x[t1] = ur + vr;
          x[t1 + 1] = ui + vi;

          x[t2] = ur - vr;
          x[t2 + 1] = ui - vi;

          t1 += mh;
          t2 = t1 + (mh << 1);

          cplxMultDiv2(&vr, &vi, x[t2 + 1], x[t2], STC(0x5a82799a),
                      STC(0x5a82799a));

          ur = x[t1] >> 1;
          ui = x[t1 + 1] >> 1;

          x[t1] = ur + vr;
          x[t1 + 1] = ui - vi;

          x[t2] = ur - vr;
          x[t2 + 1] = ui + vi;
        }
      }
      else { /* vector FFTs */
        int vl;
        const FIXP_DBL sqrt1_2 = 0x5a82799a;
        vint32m1_t vxt1p0, vxt1p1;
        vint32m1_t vxt2p0, vxt2p1;
        vint32m1_t vmul1, vmul2;
        vint32m1_t vvi, vvr, vur, vui;

        for (r = 0; r < n; r += m * vl) {
          vl = vsetvl_e32m1((n - r) >> ldm);

          // FFT round1
          vxt1p0 = vlse32_v_i32m1(&x[((r + j) << 1) + 0], 8 * m, vl);
          vxt1p1 = vlse32_v_i32m1(&x[((r + j) << 1) + 1], 8 * m, vl);
          vxt2p0 = vlse32_v_i32m1(&x[((r + j + mh) << 1) + 0], 8 * m, vl);
          vxt2p1 = vlse32_v_i32m1(&x[((r + j + mh) << 1) + 1], 8 * m, vl);

          vmul1 = vmulh_vx_i32m1(vxt2p1, sqrt1_2, vl);
          vmul2 = vmulh_vx_i32m1(vxt2p0, sqrt1_2, vl);
          vvi = vsub_vv_i32m1(vmul1, vmul2, vl);
          vvr = vadd_vv_i32m1(vmul1, vmul2, vl);

          vur = vsra_vx_i32m1(vxt1p0, 1, vl);
          vui = vsra_vx_i32m1(vxt1p1, 1, vl);

          vxt1p0 = vadd_vv_i32m1(vur, vvr, vl);
          vxt1p1 = vadd_vv_i32m1(vui, vvi, vl);

          vxt2p0 = vsub_vv_i32m1(vur, vvr, vl);
          vxt2p1 = vsub_vv_i32m1(vui, vvi, vl);

          vsse32_v_i32m1(&x[((r + j) << 1) + 0], 8 * m, vxt1p0, vl);
          vsse32_v_i32m1(&x[((r + j) << 1) + 1], 8 * m, vxt1p1, vl);
          vsse32_v_i32m1(&x[((r + j + mh) << 1) + 0], 8 * m, vxt2p0, vl);
          vsse32_v_i32m1(&x[((r + j + mh) << 1) + 1], 8 * m, vxt2p1, vl);

          // FFT round2
          vxt1p0 = vlse32_v_i32m1(&x[((r + j) << 1) + 0 + mh], 8 * m, vl);
          vxt1p1 = vlse32_v_i32m1(&x[((r + j) << 1) + 1 + mh], 8 * m, vl);
          vxt2p0 = vlse32_v_i32m1(&x[((r + j + mh) << 1) + 0 + mh], 8 * m, vl);
          vxt2p1 = vlse32_v_i32m1(&x[((r + j + mh) << 1) + 1 + mh], 8 * m, vl);

          vmul1 = vmulh_vx_i32m1(vxt2p1, sqrt1_2, vl);
          vmul2 = vmulh_vx_i32m1(vxt2p0, sqrt1_2, vl);
          vvr = vsub_vv_i32m1(vmul1, vmul2, vl);
          vvi = vadd_vv_i32m1(vmul1, vmul2, vl);

          vur = vsra_vx_i32m1(vxt1p0, 1, vl);
          vui = vsra_vx_i32m1(vxt1p1, 1, vl);

          vxt1p0 = vadd_vv_i32m1(vur, vvr, vl);
          vxt1p1 = vsub_vv_i32m1(vui, vvi, vl);

          vxt2p0 = vsub_vv_i32m1(vur, vvr, vl);
          vxt2p1 = vadd_vv_i32m1(vui, vvi, vl);

          vsse32_v_i32m1(&x[((r + j) << 1) + 0 + mh], 8 * m, vxt1p0, vl);
          vsse32_v_i32m1(&x[((r + j) << 1) + 1 + mh], 8 * m, vxt1p1, vl);
          vsse32_v_i32m1(&x[((r + j + mh) << 1) + 0 + mh], 8 * m, vxt2p0, vl);
          vsse32_v_i32m1(&x[((r + j + mh) << 1) + 1 + mh], 8 * m, vxt2p1, vl);
        }
      }
    } /* end of block 2 */
  }
}
#endif /* ifndef FUNCTION_dit_fft */

#endif /* defined(__riscv_xthead) && __riscv_v == 7000 */