//
//  MPRealFunctions.swift
//  MPFUNSwift
//
//  Created by Mike Griebling on 2022-01-24.
//

import Foundation

extension MPReal {
    
    ///   This computes the exponential function of the MPR number A and returns
    ///   the MPR result in B.  Log(2) must be precomputed to at least MPNW words
    ///   precision and the stored in the array MPLOG2CON in module MPFUNA.
    ///
    ///   This routine uses a modification of the Taylor series for Exp (t):
    ///
    ///   Exp (t) =  (1 + r + r^2 / 2! + r^3 / 3! + r^4 / 4! ...) ^ q * 2 ^ n
    ///
    ///   where the argument T has been reduced to within the closest factor of Log(2).
    ///   To further accelerate the series, the reduced argument is divided by 2^NQ.
    ///   After convergence of the series, the result is squared NQ times.  NQ is
    ///   set to bb^0.4, where bb is the number of bits of precision.
    static func mpexp (_ a:mp_real, _ b: inout mp_real, _ mpnw:Int) {
        let itrmx = 1000000, cl2 = 0.69314718055994530

        var s0 = mp_init(size: mpnw+7), s1 = mp_init(size: mpnw+7), s2 = mp_init(size: mpnw+7), s3 = mp_init(size: mpnw+7)
        var f = mp_init(size: 9)
        
        if (mpnw < 4 || a[0] < mpnw + 4 || a[0] < abs (a[2]) + 4 || b[0] < mpnw + 6) {
            print("*** MPEXP: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        var t1 = 0.0, n1 = 0
        mpmdc(a, &t1, &n1, mpnw)
        
        ///   Check for overflows and underflows.
        if n1 > 30 {
            if t1 > 0.0 {
                print("*** MPEXP: Argument is too large.")
                mpabrt (34)
            } else {
                b[1] = mpnw
                b[2] = 0
                b[3] = 0
                b[4] = 0
                b[5] = 0
                return // goto 130
            }
        }
        
        t1 = t1 * Double.pow(2.0, Double(n1))
        if abs(t1) > 1488522236.0 {
            if t1 > 0 {
                print("*** MPEXP: Argument is too large.")
                mpabrt (34)
            } else {
                b[1] = mpnw
                b[2] = 0
                b[3] = 0
                b[4] = 0
                b[5] = 0
                return // goto 130
            }
        }
        
        s0[0] = mpnw + 7
        s1[0] = mpnw + 7
        s2[0] = mpnw + 7
        s3[0] = mpnw + 7
        let mpnw1 = mpnw + 1
        
        ///   Set f1 = 1.
        f[0] = 9
        f[1] = mpnw1
        f[2] = 1
        f[3] = 0
        f[4] = 1
        f[5] = 0
        f[6] = 0
        
        ///   Check if Log(2) has been precomputed.
        var d1 = 0.0
        mpmdc(_log2, &d1, &n1, mpnw1)
        d1 = d1 * Double.pow(2.0, Double(n1))
        if (abs(d1 - cl2) > mprdfz || mpnw1 > mplog2con[1]) {
            print("*** MPEXP: Log(2) must be precomputed to precision \(mpnw1) words.",
                  "See documentation for details.")
            mpabrt (35)
        }
        
        ///   Compute the reduced argument A' = A - Log(2) * Nint [A / Log(2)].  Save
        ///   NZ = Nint [A / Log(2)] for correcting the exponent of the final result.
        mpdiv(a, _log2, &s0, mpnw1)
        mpnint(s0, &s1, mpnw1)
        mpmdc(s1, &t1, &n1, mpnw1)
        let nz = Int((t1 * Double.pow(2.0, Double(n1))).rounded())
        mpmul (_log2, s1, &s2, mpnw1)
        mpsub (a, s2, &s0, mpnw1)
        
        ///   Check if the reduced argument is zero.
        if (s0[2] == 0) {
            s0[1] = mpnw1
            s0[2] = 1
            s0[3] = 0
            s0[4] = 1
            s0[5] = 0
            s0[6] = 0
            return // goto 120
        }
        
        ///   Divide the reduced argument by 2 ^ NQ.
        let nq = max(Int(Double.pow((Double(mpnw * mpnbt)), 0.4).rounded()), 1)
        mpdivd (s0, Double.pow(2.0, Double(nq)), &s1, mpnw1)
        
        ///   Compute Exp using the usual Taylor series.
        mpeq (f, &s2, mpnw1)
        mpeq (f, &s3, mpnw1)
        var mpnw2 = mpnw1
        
        ///   The working precision used to compute each term can be linearly reduced
        ///   as the computation proceeds.
        for j in 1...itrmx {
            let t2 = Double(j)
            mpmul (s2, s1, &s0, mpnw2)
            mpdivd (s0, t2, &s2, mpnw2)
            mpadd (s3, s2, &s0, mpnw1)
            mpeq (s0, &s3, mpnw1)
            
            ///   Check for convergence of the series, and adjust working precision
            ///   for the next term.
            if (s2[2] == 0 || s2[3] < s0[3] - mpnw1) { break /*goto 100*/ }
            if j == itrmx {
                print("*** MPEXP: Iteration limit exceeded.")
                mpabrt(36)
            }
            mpnw2 = min(max(mpnw1 + (s2[3] - s0[3]) + 1, 4), mpnw1)
        }
        // 100 continue
        
        ///   Raise to the (2 ^ NQ)-th power.
        for _ in 1...nq {
            mpmul(s0, s0, &s1, mpnw1)
            mpeq(s1, &s0, mpnw1)
        }
        
        ///   Multiply by 2 ^ NZ.
        // 120 continue
        mpdmc (1.0, nz, &s2, mpnw1)
        mpmul (s0, s2, &s1, mpnw1)
        
        ///   Restore original precision level.
        mproun (&s1, mpnw)
        mpeq (s1, &b, mpnw)
        // 130 continue
    } // mpexp
    
    ///   This performs the arithmetic-geometric mean (AGM) iterations on A and B.
    ///   The AGM algorithm is as follows: Set a_0 = a and b_0 = b, then iterate
    ///
    ///    a_{k+1} = (a_k + b_k)/2
    ///    b_{k+1} = sqrt (a_k \* b_k)
    ///
    ///   until convergence (i.e., until a_k = b_k to available precision).
    ///   The result is returned in C.
    static func mpagmr (_ a:mp_real, _ b:mp_real, _ c: inout mp_real, _ mpnw:Int) {
        let itrmx = 100
        var s0 = mp_init(size: mpnw+7), s1 = mp_init(size: mpnw+7), s2 = mp_init(size: mpnw+7), s3 = mp_init(size: mpnw+7)
        
        if (mpnw < 4 || a[0] < mpnw + 4 || b[0] < abs(a[2]) + 4 || c[0] < mpnw + 6) {
            print("*** MPAGMR: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        s0[0] = mpnw + 7
        s1[0] = mpnw + 7
        s2[0] = mpnw + 7
        s3[0] = mpnw + 7
        let mpnw1 = mpnw + 1
        mpeq(a, &s1, mpnw1)
        mpeq(b, &s2, mpnw1)
        
        for j in 1...itrmx {
            mpadd (s1, s2, &s0, mpnw1)
            mpmuld (s0, 0.5, &s3, mpnw1)
            mpmul (s1, s2, &s0, mpnw1)
            mpsqrt (s0, &s2, mpnw1)
            mpeq (s3, &s1, mpnw1)
            
            ///   Check for convergence.
            mpsub (s1, s2, &s0, mpnw1)
            if (s0[2] == 0 || s0[3] < 1 - mpnw1) { break /*goto 100*/ }
            if j == itrmx {
                print("*** MPAGMR: Iteration limit exceeded.")
                mpabrt (5)
            }
        }
        // 100 continue
        mproun (&s1, mpnw)
        mpeq (s1, &c, mpnw)
    } // mpagmr
    
    ///  This computes C = A ^ B, where A, B and C are MPR.  It first checks if
    ///  B is the quotient of two integers up to 10^7 in size, in which case it
    ///  calls MPNPWR and MPNRTR.  Otherwise it calls MPLOG and MPEXP.
    static func mppower (_ a:mp_real, _ b:mp_real, _ c: inout mp_real, _ mpnw:Int) {
        let mprxx = 5.0e-10
        var s0 = mp_init(size: mpnw+7), s1 = mp_init(size: mpnw+7), s2 = mp_init(size: mpnw+7)
        
        if (mpnw < 4 || a[0] < mpnw + 4 || b[0] < abs(a[2]) + 4 || c[0] < mpnw + 6) {
            print("*** MPPOWER: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        ///  Check if A <= 0 (error), or A = 1 or B = 0 or B = 1.
        if (a[2] <= 0) {
            print("*** MPPOWER: A^B, where A is less than zero.")
            mpabrt (61)
        } else if ((a[2] == 1 && a[3] == 0 && a[4] == 1) || b[2] == 0) {
            c[1] = mpnw
            c[2] = 1
            c[3] = 0
            c[4] = 1
            c[5] = 0
            c[6] = 0
            return // goto 200
        } else if (b[2] == 1 && b[3] == 0 && b[4] == 1) {
            mpeq (a, &c, mpnw)
            return // goto 200
        }
        
        s0[0] = mpnw + 6
        s1[0] = mpnw + 6
        s2[0] = mpnw + 6
        
        ///  Check if B is rational using the extended Euclidean algorithm in DP.
        var t1 = 0.0, n1 = 0
        mpmdc (b, &t1, &n1, mpnw)
        
        var t0 = 0.0, a3 = 0.0, a4 = 0.0
        if (n1 >= -mpnbt && n1 <= mpnbt) {
            t0 = abs (t1 * Double.pow(2.0, Double(n1)))
            var t1 = max (t0, 1.0)
            var t2 = min (t0, 1.0)
            var a1 = 1.0
            var a2 = 0.0
            a3 = 0.0
            a4 = 1.0
            var exitFlag = false
            for _ in 1...20 {
                let q1 = aint (t1 / t2)
                let a5 = a1 - q1 * a3
                let a6 = a2 - q1 * a4
                let t3 = t2
                t2 = t1 - q1 * t2
                t1 = t3
                a1 = a3
                a2 = a4
                a3 = a5
                a4 = a6
                if t2 < mprxx { exitFlag = true; break /* goto 100 */ }
            }
            
            /// Call mplog and mpexp.
            if !exitFlag {
                mplog (a, &s0, mpnw)
                mpmul (s0, b, &s1, mpnw)
                mpexp (s1, &c, mpnw)
                return
            }
        }
        
        // 100 continue
        a3 = abs(a3)
        a4 = abs(a4)
        
        /// If b = a3/a4 or a4/a3 (except for sign) or then call mpnpwr and mpnrtr.
        if (abs(t0 - a3 / a4) / t0 < mprdfz) {
            a3 = sign (a3, Double(b[2]))
            mpdmc (a3, 0, &s0, mpnw)
            mpdmc (a4, 0, &s1, mpnw)
            mpdiv (s0, s1, &s2, mpnw)
            mpsub (b, s2, &s0, mpnw)
            if (s0[2] == 0 || s0[3] < b[3] + 1 - mpnw) {
                mpnpwr (a, Int(a3), &s0, mpnw)
                mpnrtr (s0, Int(a4), &c, mpnw)
                return
            }
        } else if (abs(t0 - a4 / a3) / t0 < mprdfz) {
            a4 = sign (a4, Double(b[2]))
            mpdmc (a4, 0, &s0, mpnw)
            mpdmc (a3, 0, &s1, mpnw)
            mpdiv (s0, s1, &s2, mpnw)
            mpsub (b, s2, &s0, mpnw)
            if (s0[2] == 0 || s0[3] < b[3] + 1 - mpnw) {
                mpnpwr (a, Int(a4), &s0, mpnw)
                mpnrtr (s0, Int(a3), &c, mpnw)
                return
            }
        }
        
        // 110 continue
        
        /// Call mplog and mpexp.
        mplog (a, &s0, mpnw)
        mpmul (s0, b, &s1, mpnw)
        mpexp (s1, &c, mpnw)
    } // mppower

//        subroutine mpcagm (a, b, c, mpnw)
//
//        ///   This performs the arithmetic-geometric mean (AGM) iterations on A and B
//        ///   for MPC arguments A and B.
//        ///   The AGM algorithm is as follows: Set a_0 = a and b_0 = b, then iterate
//
//        ///    a_{k+1} = (a_k + b_k)/2
//        ///    b_{k+1} = sqrt (a_k * b_k)
//
//        ///   until convergence (i.e., until a_k = b_k to available precision).
//        ///   The result is returned in C.
//
//        implicit none
//        integer, intent(in):: mpnw
//        integer itrmx, j, la, lb, lc, mp7, mpnw1
//        parameter (itrmx = 100)
//        integer (mpiknd), intent(in):: a(0:), b(0:)
//        integer (mpiknd), intent(out):: c(0:)
//        integer (mpiknd) s0(0:2*mpnw+13), s1(0:2*mpnw+13), s2(0:2*mpnw+13), &
//          s3(0:2*mpnw+13)
//
//        la = a[0]
//        lb = b[0]
//        lc = c[0]
//        if (mpnw < 4 || a[0] < abs (a[2]) + 4 || a(la) < abs (a(la+2)) + 4 &
//          || b[0] < abs (b[2]) + 4 || b(lb) < abs (b(lb+2)) + 4 || &
//          c[0] < mpnw + 6 || c(lc) < mpnw + 6) {
//          write (mpldb, 1)
//        1 format ("*** MPCAGM: uninitialized or inadequately sized arrays")
//          mpabrt (99)
//        }
//
//        mp7 = mpnw + 7
//        s0[0] = mp7
//        s0(mp7) = mp7
//        s1[0] = mp7
//        s1(mp7) = mp7
//        s2[0] = mp7
//        s2(mp7) = mp7
//        s3[0] = mp7
//        s3(mp7) = mp7
//        mpnw1 = mpnw + 1
//        mpceq (a, s1, mpnw1)
//        mpceq (b, s2, mpnw1)
//
//        for j = 1, itrmx
//          mpcadd (s1, s2, s0, mpnw1)
//          mpmuld (s0, 0.5d0, s3, mpnw1)
//          mpmuld (s0(mp7:), 0.5d0, s3(mp7:), mpnw1)
//          mpcmul (s1, s2, s0, mpnw1)
//          mpcsqrt (s0, s2, mpnw1)
//          mpceq (s3, s1, mpnw1)
//          mpcsub (s1, s2, s0, mpnw1)
//
//        ///   Check for convergence.
//
//          if ((s0[2] == 0 || s0(3) < 1 - mpnw1) &&
//            (s0(mp7+2) == 0 || s0(mp7+3) < 1 - mpnw1)) goto 100
//        }
//
//        write (mpldb, 2)
//        2 format ("*** MPCAGM: Iteration limit exceeded.")
//        mpabrt (5)
//
//        100 continue
//
//        mproun (s1, mpnw)
//        mproun (s1(mp7:), mpnw)
//        mpceq (s1, c, mpnw)
//
//        return
//        end subroutine mpcagm

    
}
