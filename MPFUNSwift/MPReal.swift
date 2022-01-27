//
//  MPReal.swift
//  MPFUNSwift
//
//  Created by Mike Griebling on 2022-01-20.
//

import Foundation
import Numerics    // Give access to a Real & Complex types

public typealias Complex = ComplexModule.Complex<Double>

/// This is a Swift port of the original author"s MPFUN20 fortran code. Why for this? Unless you"re a fortran junkie, the
/// original code is hugely unreadable (in my opinion) with no disrespect to the original author -- one can only do so much
/// in making fortran code readable.  The Swift code is, hopefully, more understandable and perhaps faster than the original
/// fortran code.  Benchmarks will follow.
///
/// From the original author"s comments:
///
/// This package permits one to perform floating-point computations (real and complex) to arbitrarily high numeric precision,
/// by making only relatively minor changes to existing Fortran-90 programs.  All basic arithmetic operations and transcendental
/// functions are supported, together with several special functions.
///
/// In addition to fast execution times, one key feature of this package is a 100% THREAD-SAFE design, which means
/// that user-level applications can be easily converted for parallel execution, say using a threaded parallel environment
/// such as OpenMP.
///
public struct MPReal : Codable, CustomStringConvertible, ExpressibleByArrayLiteral {

    // used internally for convenience
    typealias mp_real = [Int]
    
    /// The star of the show: an (n+6)-long vector of 64-bit integers, where *n* is the number of mantissa words.
    /// Words definitions:
    /// - Word 0: Total space allocated for this array in 64-bit words
    /// - Word 1: The working precision level (in bits) associated with this data.
    /// - Word 2: The sign of the value.
    /// - Word 3: The exponent, base 2.
    /// - Word 4: The first word of the mantissa
    /// - Word 5 to N + 4: Mantissa words (signed integers between 0 and 2**60 − 1).
    /// - Word N + 5: Not used at present.
    var n : [Int]
    
    // MARK: - Initialization methods

    static func mp_init(size: Int) -> [Int] { var x = [Int](repeating: 0, count: size); x[0] = size; return x }
    static func mp_zero(_ a: inout mp_real, _ mpnw:Int) { a[1...5] = [mpnw, 0, 0, 0, 0] }

    public init(_ ni: Int = 0, size: Int = 0) {
        let size = max(size, 4) // need at least 4 words
        n = MPReal.mp_init(size: size+6)  // allocates size+6 words
        n[0] = n.count           // total space
        n[1] = size              // space used
        n[2] = ni.signum()       // sign
        n[3] = 0                 // exponent
        n[4] = Int(ni.magnitude) >> MPReal.mpnbt // value
        n[5] = 0
    }
    
    /// Use only if knowledgeable about the internal structure
    public init(arrayLiteral elements: Int...) { self.n = elements }
    init(_ mp: mp_real) { self.n = mp }
    
    public init(_ d:Double) {
        self.init(0, size: 4)
        MPReal.mpdmc(d, 0, &n, 4)
    }
    
    public init(_ s:String) {
        self.init(0, size: s.count/Int(MPReal.mpdpw)+6)
        MPReal.mpctomp(s, &n, n.count-6)
    }
    
    // MARK: - Custom string convertible compliance
    public var description: String {
        let digits = Int(MPReal.mpdpw * Double(self.n[2]) + 5)
        return toString(with: digits)
    }
    
    public func toString(with digits:Int) -> String {
        var str = ""
        MPReal.mpfformat(self.n, digits+20, digits, &str, 10)
        while str.hasSuffix("0") && str.count > 1 { str.removeLast() } // remove trailing zeros
        if str.hasSuffix(".") { str.removeLast() } // remove trailing "."
        return str
    }
    
    // need to be initialized by mpinifft()
    static var mpuu1 = [Complex]()
    static var mpuu2 = [Complex]()
}

extension MPReal : ExpressibleByIntegerLiteral {
    public init(integerLiteral value: Int) { self.init(value) }
}

extension MPReal : Real {
    
    public static var pi: MPReal { MPReal(_pi) }
    
    public mutating func round(_ rule: FloatingPointRoundingRule) {
        // STUB
    }
    
    public init(sign: FloatingPointSign, exponent: Int, significand: MPReal) {
        let words = significand.n[1]
        var res = MPReal.mp_init(size: significand.n.count)
        MPReal.mpdmc(2.0, exponent, &res, words)
        MPReal.mpmul(res, significand.magnitude.n, &res, words)
        if sign == .minus { MPReal.mpneg(res, &res, words) }
        self.init(res)
    }
    
    public init(signOf x: MPReal, magnitudeOf y: MPReal) {
        let sign = x.sign
        let mag = sign == .minus ? -y.magnitude : y.magnitude
        self.init(mag.n)
    }
    
    public init<Source>(_ value: Source) where Source : BinaryInteger {
        // STUB
        self.init(Int(value), size: 0)
    }
    
    public static var radix: Int { mpbdx }
    public static var nan: MPReal { zero } // STUB
    public static var signalingNaN: MPReal { zero } // STUB
    public static var infinity: MPReal {  zero } // STUB
    
    public static var greatestFiniteMagnitude: MPReal { zero } // STUB
    
    public var ulp: MPReal { MPReal.zero } // STUB
    public static var leastNormalMagnitude: MPReal { zero } // STUB
    public static var leastNonzeroMagnitude: MPReal { zero } // STUB
    public var sign: FloatingPointSign { self.n[2] < 0 ? .minus : .plus }
    public var exponent: Int { self.n[3] }
    public var significand: MPReal { MPReal.zero } // STUB
    public mutating func formRemainder(dividingBy other: MPReal) {  } // STUB
    public mutating func formTruncatingRemainder(dividingBy other: MPReal) { } // STUB
    
    public mutating func formSquareRoot() {
        var z = MPReal.mp_init(size: self.n.count)
        MPReal.mpsqrt(self.n, &z, self.n.count-4)
        self.n = z
    }
    
    public mutating func addProduct(_ lhs: MPReal, _ rhs: MPReal) { self += lhs * rhs }
    public var nextUp: MPReal { MPReal.zero } // STUB
    
    public func isEqual(to other: MPReal) -> Bool {
        var ic = 0
        MPReal.mpcpr(self.n, other.n, &ic, self.n.count-4)
        return ic == 0
    }
    
    public func isLess(than other: MPReal) -> Bool {
        var ic = 0
        MPReal.mpcpr(self.n, other.n, &ic, self.n.count-4)
        return ic == -1
    }
    
    public func isLessThanOrEqualTo(_ other: MPReal) -> Bool {
        var ic = 0
        MPReal.mpcpr(self.n, other.n, &ic, self.n.count-4)
        return ic <= 0
    }
    
    public func isTotallyOrdered(belowOrEqualTo other: MPReal) -> Bool { self <= other }
    
    public var isNormal: Bool { true } // STUB
    public var isFinite: Bool { true } // STUB
    public var isZero: Bool { self.n[2] == 0 }
    public var isSubnormal: Bool { false } // STUB
    public var isInfinite: Bool { false } // STUB
    public var isNaN: Bool { false } // STUB
    public var isSignalingNaN: Bool { false } // STUB
    public var isCanonical: Bool { false } // STUB
    public func distance(to other: MPReal) -> MPReal { MPReal.zero } // STUB
    public func advanced(by n: MPReal) -> MPReal { self + n }
    
    public static func - (lhs: MPReal, rhs: MPReal) -> MPReal {
        var z = mp_init(size: max(lhs.n.count, rhs.n.count))
        MPReal.mpsub(lhs.n, rhs.n, &z, z.count-6)
        return MPReal(z)
    }
    
    public static func atan2(y: MPReal, x: MPReal) -> MPReal { y } // STUB
    public static func erf(_ x: MPReal) -> MPReal { x } // STUB
    public static func erfc(_ x: MPReal) -> MPReal { x } // STUB
    public static func exp2(_ x: MPReal) -> MPReal { pow(MPReal(2), x) }
    public static func hypot(_ x: MPReal, _ y: MPReal) -> MPReal { x } // STUB
    public static func gamma(_ x: MPReal) -> MPReal { x } // STUB
    public static func log2(_ x: MPReal) -> MPReal { x } // STUB
    public static func log10(_ x: MPReal) -> MPReal { x } // STUB
    public static func logGamma(_ x: MPReal) -> MPReal { x } // STUB
    
    public static func /= (a: inout MPReal, b: MPReal) {
        var z = mp_init(size: a.n.count)
        MPReal.mpdiv(a.n, b.n, &z, a.n.count-4)
        a = MPReal(z)
    }
    
    public static func exp(_ x: MPReal) -> MPReal {
        var y = mp_init(size: x.n.count)
        mpexp(x.n, &y, x.n.count-4)
        return MPReal(y)
    }
    
    public static func expMinusOne(_ x: MPReal) -> MPReal { x } // STUB
    public static func cosh(_ x: MPReal) -> MPReal { x } // STUB
    public static func sinh(_ x: MPReal) -> MPReal { x } // STUB
    public static func tanh(_ x: MPReal) -> MPReal { x } // STUB
    public static func cos(_ x: MPReal) -> MPReal { x } // STUB
    public static func sin(_ x: MPReal) -> MPReal { x } // STUB
    public static func tan(_ x: MPReal) -> MPReal { x } // STUB
    
    public static func log(_ x: MPReal) -> MPReal {
        var y = mp_init(size: x.n.count)
        mplog(x.n, &y, x.n.count-4)
        return MPReal(y)
    }
    
    public static func log(onePlus x: MPReal) -> MPReal {  x } // STUB
    public static func acosh(_ x: MPReal) -> MPReal { x } // STUB
    public static func asinh(_ x: MPReal) -> MPReal { x } // STUB
    public static func atanh(_ x: MPReal) -> MPReal { x } // STUB
    public static func acos(_ x: MPReal) -> MPReal { x } // STUB
    public static func asin(_ x: MPReal) -> MPReal { x } // STUB
    public static func atan(_ x: MPReal) -> MPReal { x } // STUB
    
    public static func pow(_ x: MPReal, _ y: MPReal) -> MPReal {
        var z = mp_init(size: x.n.count)
        mppower(x.n, y.n, &z, x.n.count-4)
        return MPReal(z)
    }
    
    public static func pow(_ x: MPReal, _ n: Int) -> MPReal {
        var z = mp_init(size: x.n.count+2)
        mpnpwr(x.n, n, &z, x.n.count-4)
        return MPReal(z)
    }
    
    public static func root(_ x: MPReal, _ n: Int) -> MPReal {
        var z = mp_init(size: x.n.count+2)
        mpnrtr(x.n, n, &z, x.n.count-4)
        return MPReal(z)
    }
    
    public init?<T>(exactly source: T) where T : BinaryInteger {
        // STUB
        self.init(Int(source), size: 0)
    }
    
    public var magnitude: MPReal {
        var z = MPReal.mp_init(size: self.n.count)
        MPReal.mpabs(self.n, &z, self.n.count-4)
        return MPReal(z)
    }
    
    public static func * (lhs: MPReal, rhs: MPReal) -> MPReal {
        var z = mp_init(size: lhs.n.count)
        MPReal.mpmul(lhs.n, rhs.n, &z, lhs.n.count-4)
        return MPReal(z)
    }
    
    public static func *= (lhs: inout MPReal, rhs: MPReal) { lhs = lhs * rhs }
    
    public static func + (lhs: MPReal, rhs: MPReal) -> MPReal {
        var z = mp_init(size: lhs.n.count)
        MPReal.mpadd(lhs.n, rhs.n, &z, lhs.n.count-4)
        return MPReal(z)
    }
    
    public static func ** (base: MPReal, power: Int) -> MPReal { pow(base, power) }
    public static func ** (base: MPReal, power: MPReal) -> MPReal { pow(base, power) }
}

// Definition of power operator
infix operator ** : ExponentPrecedence
precedencegroup ExponentPrecedence {
    associativity: left
    higherThan: MultiplicationPrecedence
}

// Constant calculations and stored values
extension MPReal {
    
    // Constants to the current precision set by mpinitran()
    static var _pi = mp_real()
    static var _log2 = mp_real()
    static var _egamma = mp_real()
    
    public static func pi(_ mpnw:Int) -> MPReal {
        if _pi.count == 0 { mpinitran(mpnw) }
        return MPReal(_pi)
    }
    public static func log2(_ mpnw:Int) -> MPReal {
        if _log2.count == 0 { mpinitran(mpnw) }
        return MPReal(_log2)
    }
    public static func egamma(_ mpnw:Int) -> MPReal {
        if _egamma.count == 0 { mpinitran(mpnw) }
        return MPReal(_egamma)
    }
    
    /// This routine computes pi, log(2) and egamma, and stores this data in the
    /// proper arrays in module MPFUNA.  MPNW is the largest precision level
    /// (in words) that will be subsequently required for this run at the user level.
    static func mpinitran(_ mpnw:Int) {
        //   Add three words to mpnw, since many of the routines in this module
        //   increase the working precision level by one word upon entry.
        let nwds = mpnw + 3
        
        //  Compute pi, log(2) and sqrt(2)/2.
        let nwds6 = nwds + 6
        _log2   = mp_init(size: nwds6)  // empty
        _pi     = mp_init(size: nwds6)  // empty
        _egamma = mp_init(size: nwds6)  // empty
        
        mppiq(&_pi, nwds)
        mplog2q(&_log2, nwds)
        mpegammaq(&_egamma, nwds - 3)
    } // mpinitran
    
    ///   This computes Euler's gamma to available precision (MPNW mantissa words).
    ///   The algorithm is the following, due to Brent and MacMillan:

    ///   Given a desired precision D digits, select n1 = 0.6*D and n2 = 2.2*D
    ///   (approximately). Set A_0 = -log(n1); B_0 = 1; U_0 = A_0; V_0 = 1.
    ///   Then iterate for k = 1, 2, ..., n2:
    ///     B_k = B_{k-1} * n1^2 / k^2
    ///     A_k = (A_{k-1} * n1^2/k + B_k) / k
    ///     U_k = U_{k-1} + A_k
    ///     V_k = V_{k-1} + B_k
    ///   Then set gamma = U_k / V_k.
    static func mpegammaq (_ egamma: inout mp_real, _ mpnw:Int) {
        var aa = mp_init(size: mpnw+7), bb = mp_init(size: mpnw+7), uu = mp_init(size: mpnw+7), vv = mp_init(size: mpnw+7)
        var s0 = mp_init(size: mpnw+7), s1 = mp_init(size: mpnw+7)
        
        if (mpnw < 4 || egamma[0] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        let mpnw1 = mpnw + 1
        let ndp = Double(mpnw - 1) * mpdpw
        let n1 = Int((0.6 * ndp).rounded())
        let n2 = Int((2.2 * ndp).rounded())
        
        ///   Unless precision is very high, just copy egamma from table.
        if mpnw1 <= mpegammacon[1] {
            mpeq(mpegammacon, &egamma, mpnw)
            return // goto 110
        }
        
        ///   Initialize aa, bb, uu, vv, s0, s1.
        aa[0] = mpnw1 + 6
        bb[0] = mpnw1 + 6
        uu[0] = mpnw1 + 6
        vv[0] = mpnw1 + 6
        s0[0] = mpnw1 + 6
        s1[0] = mpnw1 + 6
        
        /// aa = -log (mpreal (dble (n1), nwds))
        mpdmc(Double(n1), 0, &s0, mpnw1)
        mplog(s0, &aa, mpnw1)
        aa[2] = -aa[2]
        
        /// bb = mpreal (1.d0, nwds)
        /// uu = aa
        /// vv = bb
        mpdmc(1.0, 0, &bb, mpnw1)
        mpeq(aa, &uu, mpnw1)
        mpeq(bb, &vv, mpnw1)
        
        for k in 1...n2 {
            ///   bb = bb * dble (n1)**2 / dble (k)**2
            let kd = Double(k)
            let x2 = Double(n1) * Double(n1)
            let k2 = kd * kd
            mpmuld(bb, x2, &s0, mpnw1)
            mpdivd(s0, k2, &bb, mpnw1)
            
            ///  aa = (aa * dble (n1)**2 / dble (k) + bb) / dble (k)
            mpmuld(aa, x2, &s0, mpnw1)
            mpdivd(s0, kd, &s1, mpnw1)
            mpadd(s1, bb, &s0, mpnw1)
            mpdivd(s0, kd, &aa, mpnw1)
            
            ///  uu = uu + aa
            ///  vv = vv + bb
            mpadd(uu, aa, &s0, mpnw1)
            mpeq(s0, &uu, mpnw1)
            mpadd(vv, bb, &s0, mpnw1)
            mpeq(s0, &vv, mpnw1)
            
            ///  if (aa%mpr(3) < uu%mpr(3) - nwds && bb%mpr(3) < vv%mpr(3) - nwds) &
            if (aa[3] < uu[3] - mpnw1 && bb[3] < vv[3] - mpnw1) { break /* goto 100 */ }
        }
        
        // 100 continue
        
        /// egammax = uu / vv
        mpdiv(uu, vv, &s0, mpnw1)
        
        ///   Restore original precision level.
        mproun(&s0, mpnw)
        mpeq(s0, &egamma, mpnw)
    } // mpegammaq
    
    ///   This computes the natural logarithm of the MPR number A and returns the MPR
    ///   result in B.
    
    ///   The Taylor series for Log converges much more slowly than that of Exp.
    ///   Thus this routine does not employ Taylor series (except if the argument
    ///   is extremely close to 1), but instead computes logarithms by solving
    ///   Exp (b) = a using the following Newton iteration:
    
    ///     x_{k+1} = x_k + [a - Exp (x_k)] / Exp (x_k)
    
    ///   These iterations are performed with a maximum precision level MPNW that
    ///   is dynamically changed, approximately doubling with each iteration.
    static func mplog (_ a: mp_real, _ b: inout mp_real, _ mpnw:Int) {
        let alt = 0.693147180559945309, cl2 = 1.4426950408889633
        let rtol = Foundation.pow(0.5, Double(7)), itrmax = 1000000, nit = 3
 
        var s0 = mp_init(size: mpnw+7), s1 = mp_init(size: mpnw+7), s2 = mp_init(size: mpnw+7), s3 = mp_init(size: mpnw+7)
        var s4 = mp_init(size: mpnw+7), f1 = mp_init(size: 9)
        
        if (mpnw < 4 || a[0] < mpnw + 4 || a[0] < abs(a[2]) + 4 || b[0] < mpnw + 6) {
            print("*** MPLOG: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        let ia = sign(1, a[2])
        let na = min(abs(a[2]), mpnw)
        
        if (ia < 0 || na == 0) {
            print("*** MPLOG: Argument is less than or equal to zero.")
            mpabrt (50)
        }
        
        ///   Check if input is exactly one.
        if (a[2] == 1 && a[3] == 0 && a[4] == 1) {
//            b[1] = mpnw
//            b[2] = 0
//            b[3] = 0
//            b[4] = 0
//            b[5] = 0
            mp_zero(&b, mpnw)
            return // goto 130
        }
        
        var mpnw1 = mpnw + 1
        s0[0] = mpnw + 7
        s1[0] = mpnw + 7
        s2[0] = mpnw + 7
        s3[0] = mpnw + 7
        s4[0] = mpnw + 7
        
        f1[0] = 9
        f1[1] = mpnw1
        f1[2] = 1
        f1[3] = 0
        f1[4] = 1
        f1[5] = 0
        f1[6] = 0
        
        ///   If the argument is sufficiently close to 1, employ a Taylor series.
        mpsub (a, f1, &s0, mpnw1)
        if (s0[2] == 0 || s0[3] <= min(-2, Int(-rtol * Double(mpnw1)))) {
            mpeq(s0, &s1, mpnw1)
            mpeq(s1, &s2, mpnw1)
            // let i1 = 1
            var iss = 1
            let tol = s0[3] - mpnw1
            
            for i1 in 2...itrmax {
                iss = -iss
                let st = Double(iss * i1)
                mpmul(s1, s2, &s3, mpnw1)
                mpeq(s3, &s2, mpnw1)
                mpdivd(s3, st, &s4, mpnw1)
                mpadd(s0, s4, &s3, mpnw1)
                mpeq(s3, &s0, mpnw1)
                if (s4[2] == 0 || s4[3] < tol) {
                    mproun(&s3, mpnw)
                    mpeq(s3, &b, mpnw)
                    return // goto 120
                }
            }
            print("*** MPLOG: Iteration limit exceeded \(itrmax)")
            mpabrt (54)
        }
        
        ///   Determine the least integer MQ such that 2 ^ MQ .GE. MPNW.
        let t2 = Double(mpnw)
        let mq = Int(cl2 * Foundation.log(t2) + 2.0 - mprdfz)
        
        ///   Compute initial approximation of Log (A).
        var n1 = 0, t1 = 0.0
        mpmdc(a, &t1, &n1, mpnw)
        t1 = Foundation.log(t1) + Double(n1) * alt
        mpdmc (t1, 0, &s3, mpnw)
        mpnw1 = 4
        var iq = 0
        
        ///   Perform the Newton-Raphson iteration described above with a dynamically
        ///   changing precision level MPNW (one greater than powers of two).
        for k in 0...mq {
            if (k > 1) { mpnw1 = min(2 * mpnw1 - 2, mpnw) + 1 }
            
            // 110  continue
            while true {
                mpexp (s3, &s0, mpnw1)
                mpsub (a, s0, &s1, mpnw1)
                mpdiv (s1, s0, &s2, mpnw1)
                mpadd (s3, s2, &s1, mpnw1)
                mpeq (s1, &s3, mpnw1)
                if (k == mq - nit && iq == 0) {
                    iq = 1
//                    goto 110
                }
                break
            }
        }
        
        ///   Restore original precision level.
        mproun (&s3, mpnw)
        mpeq (s3, &b, mpnw)
        
        // 130 continue
    } // mplog

    ///   This computes log(2) to mpnw words precision, using an algorithm due to Salamin
    ///   and Brent:  Select n > 2^m, where m is the number of bits of desired precision
    ///   precision in the result.  Then
    ///
    ///   Log(2) = Pi / [2 AGM (1, 4/x)]
    ///
    ///   Where AGM (a, b) denotes the arithmetic-geometric mean:  Set a_0 = a and
    ///   b_0 = b, then iterate
    ///    a_{k+1} = (a_k + b_k)/2
    ///    b_{k+1} = sqrt (a_k * b_k)
    ///   until convergence (i.e., until a_k = b_k to available precision).
    static func mplog2q (_ alog2: inout mp_real, _ mpnw:Int) {
        let cpi = 3.141592653589793238
        var s1 = mp_init(size: mpnw+7), s2 = mp_init(size: mpnw+7), s3 = mp_init(size: mpnw+7), s4 = mp_init(size: mpnw+7)
        var f1 = mp_init(size: 9), f4 = mp_init(size: 9)
        
        /// End of declaration
        
        if (mpnw < 4 || alog2[0] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        ///   Define sections of the scratch array.
        s1[0] = mpnw + 7
        s2[0] = mpnw + 7
        s3[0] = mpnw + 7
        s4[0] = mpnw + 7
        var mpnw1 = mpnw + 1
        
        ///   Unless precision is very high, just copy log2 from table.
        if mpnw1 <= mplog2con[1] {
            mpeq(mplog2con, &alog2, mpnw)
            return // goto 100
        }
        
        ///   Check if Pi has been precomputed.
        var d1 = 0.0; var n1 = 0
        mpmdc(_pi, &d1, &n1, mpnw)
        d1 = d1 * Double(1<<n1) // 2.0**n1
        if (abs(d1 - cpi) > mprdfz || mppicon[1] < mpnw) {
            print("*** \(#function): Pi must be precomputed to precision \(mpnw) words.",
                      "See documentation for details.")
            mpabrt (53)
        }
        
        ///   Define sections of the scratch array.
        s1[0] = mpnw + 7
        s2[0] = mpnw + 7
        s3[0] = mpnw + 7
        s4[0] = mpnw + 7
        mpnw1 = mpnw + 1
        
        ///   Set f1 = 1.
        f1[0] = 9
        f1[1] = mpnw1
        f1[2] = 1
        f1[3] = 0
        f1[4] = 1
        f1[5] = 0
        f1[6] = 0
        
        ///   Set f4 = 4.
        f4[0] = 9
        f4[1] = mpnw1
        f4[2] = 1
        f4[3] = 0
        f4[4] = 4
        f4[5] = 0
        f4[6] = 0
        
        ///   Set s4 to 2^(n/2), where n is the number of bits desired. n48 = n/mpnbt.
        ///   Note that this value can be directly set in the first few words of s4,
        ///   avoiding explicit exponentiation.
        let n = mpnbt * (mpnw1 / 2 + 2)
        let n48 = n / mpnbt
        
        s4[1] = mpnw1
        s4[2] = 1
        s4[3] = n48
        s4[4] = 1
        s4[5] = 0
        s4[6] = 0
        
        ///   Perform AGM iterations.
        mpeq (f1, &s1, mpnw1)
        mpdiv (f4, s4, &s2, mpnw1)
        mpagmr (s1, s2, &s3, mpnw1)
        
        ///   Set Log(2) = Pi / (2 * N * S3), where S3 is the limit of the AGM iterations.
        mpmuld (s3, 2.0 * Double(n), &s1, mpnw1)
        mpdiv (_pi, s1, &s2, mpnw1)
        mproun (&s2, mpnw)
        mpeq (s2, &alog2, mpnw)
        //    100 continue
    } // mplog2q
        
    ///   This computes Pi to available precision (MPNW mantissa words).
    ///   The algorithm that is used for computing Pi, which is due to Salamin
    ///   and Brent, is as follows:
    ///
    ///   Set  A_0 = 1,  B_0 = 1/Sqrt(2)  and  D_0 = Sqrt(2) - 1/2.
    ///
    ///   Then from k = 1 iterate the following operations:
    ///
    ///   A_k = 0.5 * (A_{k-1} + B_{k-1})
    ///   B_k = Sqrt (A_{k-1} * B_{k-1})
    ///   D_k = D_{k-1} - 2^k * (A_k - B_k) ^ 2
    ///
    ///   Then  P_k = (A_k + B_k) ^ 2 / D_k  converges quadratically to Pi.
    ///   In other words, each iteration approximately doubles the number of correct
    ///   digits, providing all iterations are done with the maximum precision.
    ///   The constant cl2 (below) = 1 / log(2) (DP approximation).
    static func mppiq (_ pi: inout mp_real, _ mpnw:Int) {
        let cl2 = 1.4426950408889633  // 1 / log(2)
        var s0 = mp_init(size: mpnw+7), s1 = mp_init(size: mpnw+7), s2 = mp_init(size: mpnw+7)
        var s3 = mp_init(size: mpnw+7), s4 = mp_init(size: mpnw+7), f = mp_init(size: 9)
        
        if mpnw < 4 || pi[0] < mpnw + 6 {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        s0[0] = mpnw + 7
        s1[0] = mpnw + 7
        s2[0] = mpnw + 7
        s3[0] = mpnw + 7
        s4[0] = mpnw + 7
        let mpnw1 = mpnw + 1
        
        ///   Unless precision is very high, just copy pi from table.
        if mpnw1 <= mppicon[1] {
            mpeq(mppicon, &pi, mpnw)
            return
            // goto 100
        }
        
        ///   Determine the number of iterations required for the given precision level.
        ///   This formula is good only for this Pi algorithm.
        var t1 = Double(mpnw1) * Double.log10(Double(mpbdx))
        let mq = Int(cl2 * (Double.log(t1) - 1.0) + 1.0)
        
        ///   Initialize as above.
        s0[1] = mpnw
        s0[2] = 1
        s0[3] = 0
        s0[4] = 1
        s0[5] = 0
        s0[6] = 0
        f[0] = 9
        f[1] = mpnw1
        f[2] = 1
        f[3] = 0
        f[4] = 2
        f[5] = 0
        f[6] = 0
        mpsqrt(f, &s2, mpnw1)
        mpmuld(s2, 0.5, &s1, mpnw1)
        f[3] = -1
        f[4] = Int(0.5 * Double(mpbdx))
        mpsub(s2, f, &s4, mpnw1)
        
        ///   Perform iterations as described above.
        for k in 1...mq {
            mpadd(s0, s1, &s2, mpnw1)
            mpmul(s0, s1, &s3, mpnw1)
            mpsqrt(s3, &s1, mpnw1)
            mpmuld(s2, 0.5, &s0, mpnw1)
            mpsub(s0, s1, &s2, mpnw1)
            mpmul(s2, s2, &s3, mpnw1)
            t1 = Double.pow(2.0, Double(k))
            mpmuld(s3, t1, &s2, mpnw1)
            mpsub(s4, s2, &s3, mpnw1)
            mpeq(s3, &s4, mpnw1)
        }
        
        ///   Complete computation.
        mpadd(s0, s1, &s2, mpnw1)
        mpmul(s2, s2, &s2, mpnw1)
        mpdiv(s2, s4, &s2, mpnw1)
        mpeq(s2, &s0, mpnw1)
        
        ///   Restore original precision level.
        mproun(&s0, mpnw)
        mpeq(s0, &pi, mpnw)
        
        // 100 continue
    } // mppiq
    
}

/// Primarily converts to/from strings and allows I/O
extension MPReal {
    
    // MARK: - Elementary operations
    
    static func mpabrt( _ ier:Int) { assertionFailure("***  \(#function): Execution terminated, error code = \(ier)") }
    
    ///   This routine sets rb = absolute value of ra.
    static func mpabs (_ ra:mp_real, _ rb: inout mp_real, _ mpnw:Int) {
        mpeq (ra, &rb, mpnw)  // rb = ra
        rb[2] = min(Int(abs(ra[2])), mpnw)
    } // mpabs
    
    ///   This routine adds MPR numbers A and B to yield C.
    static func mpadd (_ a: mp_real, _ b: mp_real, _ c: inout mp_real, _ mpnw:Int) {
        
        if (mpnw < 4 || a[0] < abs(a[2]) + 4 || b[0] < abs(b[2]) + 4 || c[0] < mpnw + 6) {
            print("***  \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        var d = mp_init(size: mpnw+7)
        
        let ia = sign(1, a[2])
        let ib = sign(1, b[2])
        let na = min(abs(a[2]), mpnw)
        let nb = min(abs(b[2]), mpnw)
        
        ///   Check for zero inputs.
        if na == 0 {
            ///   A is zero -- the result is B.
            c[1] = mpnw
            c[2] = sign(nb, ib)
            c[3...nb+4] = b[3...nb+4]
//            for i in 2...nb+3 {
//                c[i+1] = b[i+1]
//            }
            c[nb+4] = 0
            c[nb+5] = 0
            return
        } else if nb == 0 {
            ///   B is zero -- the result is A.
            c[1] = mpnw
            c[2] = sign(na, ia)
            c[3...na+4] = a[3...na+4]
//            for i in 2...na+3 {
//                c[i+1] = a[i+1]
//            }
//
            c[na+4] = 0
            c[na+5] = 0
            return
        }
        
        let idb:Int
        if ia == ib {
            idb = 1
        } else {
            idb = -1
        }
        let ixa = a[3]
        let ixb = b[3]
        let ish = ixa - ixb
        let ixd, nd:Int
        if ish >= 0 {
            ///   A has greater exponent than B, so B must be shifted to the right.
            ///  m1 = number of initial A words to be copied without adding B.
            ///  m2 = end index of A words to be added to shifted B words, after end of initial A.
            ///  m3 = end index of A words to be copied without adding, after end of shifted B section.
            ///  m4 = end index of zero words after the end of A words.
            ///  m5 = end index of B words to be copied with a shift, after end of A words.
            let m1 = min (na, ish)
            let m2 = min (na, nb + ish)
            let m3 = na
            let m4 = min (max (na, ish), mpnw + 1)
            let m5 = min (max (na, nb + ish), mpnw + 1)
            
            for i in stride(from: 1, through: m1, by:1) { //1...m1 {
                d[i+3] = a[i+3]
            }
            
            for i in stride(from: m1+1, through: m2, by: 1) {
                d[i+3] = a[i+3] + idb * b[i+2-ish+1]
            }
            
            for i in stride(from: m2+1, through: m3, by: 1) {
                d[i+3] = a[i+3]
            }
            
            for i in stride(from: m3+1, through: m4, by: 1) {
                d[i+3] = 0
            }
            
            for i in stride(from: m4+1, through: m5, by: 1) {
                d[i+3] = idb * b[i+2-ish+1]
            }
            
            nd = m5
            ixd = ixa
            d[nd+4] = 0
            d[nd+5] = 0
        } else {
            ///   B has greater exponent than A, so A must be shifted to the right.
            let nsh = -ish
            let m1 = min(nb, nsh)
            let m2 = min(nb, na + nsh)
            let m3 = nb
            let m4 = min(max(nb, nsh), mpnw + 1)
            let m5 = min(max(nb, na + nsh), mpnw + 1)
            
            for i in 1...m1 {
                d[i+3] = idb * b[i+3]
            }
            
            for i in stride(from: m1+1, through: m2, by: 1) {
                d[i+3] = a[i+2-nsh+1] + idb * b[i+3]
            }
            
            for i in stride(from: m2+1, through: m3, by: 1) {
                d[i+3] = idb * b[i+3]
            }
            
            for i in stride(from: m3+1, through: m4, by: 1) {
                d[i+3] = 0
            }
            
            for i in stride(from: m4+1, through: m5, by: 1) {
                d[i+3] = a[i+2-nsh+1]
            }
            
            nd = m5
            ixd = ixb
            d[nd+4] = 0
            d[nd+5] = 0
        }
        
        ///   Call mpnorm to fix up result and store in c.
        d[0] = mpnw + 6
        d[1] = mpnw
        d[2] = sign(nd, ia)
        d[3] = ixd
        mpnorm(&d, &c, mpnw)
    } // mpadd

    // MARK: - Complex operations
    ///   This routine returns the absolute value of the MPC argument A (the
    ///   result is of type MPR).
    static func mpcabs (_ a:mp_real, _ b: inout mp_real, _ mpnw:Int) {
        let n = mpnw+7
        var s0 = mp_init(size: n), s1 = mp_init(size: n), s2 = mp_init(size: n)
        
        let la = a[0]
        if (mpnw < 4 || a[0] < abs(a[2]) + 4 || a[la] < abs(a[la+2]) + 4 || b[0] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        let mpnw1 = mpnw + 1
        s0[0] = mpnw + 7
        s1[0] = mpnw + 7
        s2[0] = mpnw + 7
        mpmul(a, a, &s0, mpnw1)
        mpmul(mp_real(a[la...]), mp_real(a[la...]), &s1, mpnw1)
        mpadd(s0, s1, &s2, mpnw1)
        mpsqrt(s2, &s0, mpnw1)
        mproun(&s0, mpnw)
        mpeq(s0, &b, mpnw)
    } // mpcabs

    ///   This routine adds the MPC numbers A and B.
    static func mpcadd (_ a:mp_real, _ b:mp_real, _ c:inout mp_real, _ mpnw:Int) {
        let la = a[0]
        let lb = b[0]
        let lc = c[0]
        if (mpnw < 4 || a[0] < abs(a[2]) + 4 || a[la] < abs(a[la+2]) + 4 ||
            b[0] < abs(b[2]) + 4 || b[lb] < abs(b[lb+2]) + 4 ||
            c[0] < mpnw + 6 || c[lc] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        mpadd(a, b, &c, mpnw)
        var t = mp_real(c[lc...])
        mpadd(mp_real(a[la...]), mp_real(b[lb...]), &t, mpnw)
        c[lc...] = t[0...]
    } // mpcadd

    ///   This routine divides the MPC numbers A and B.
    static func mpcdiv (_ a:mp_real, _ b:mp_real, _ c:inout mp_real, _ mpnw:Int) {
        var s0 = mp_init(size: mpnw+7), s1 = mp_init(size: mpnw+7), s2 = mp_init(size: mpnw+7)
        var s3 = mp_init(size: mpnw+7), s4 = mp_init(size: mpnw+7)
        let la = a[0]
        let lb = b[0]
        let lc = c[0]
        if (mpnw < 4 || a[0] < abs(a[2]) + 4 || a[la] < abs(a[la+2]) + 4 ||
            b[0] < abs(b[2]) + 4 || b[lb] < abs(b[lb+2]) + 4 ||
            c[0] < mpnw + 6 || c[lc] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        let mpnw1 = mpnw + 1
        s0[0] = mpnw + 7
        s1[0] = mpnw + 7
        s2[0] = mpnw + 7
        s3[0] = mpnw + 7
        s4[0] = mpnw + 7
        
        mpmul(a, b, &s0, mpnw1)
        mpmul(mp_real(a[la...]), mp_real(b[lb...]), &s1, mpnw1)
        mpadd(s0, s1, &s2, mpnw1)
        mpmul(a, mp_real(b[lb...]), &s0, mpnw1)
        s0[2] = -s0[2]
        mpmul(mp_real(a[la...]), b, &s1, mpnw1)
        mpadd(s0, s1, &s3, mpnw1)
        
        mpmul(b, b, &s0, mpnw1)
        mpmul(mp_real(b[lb...]), mp_real(b[lb...]), &s1, mpnw1)
        mpadd(s0, s1, &s4, mpnw1)
        mpdiv(s2, s4, &s0, mpnw1)
        mpdiv(s3, s4, &s1, mpnw1)
        
        mproun(&s0, mpnw)
        mproun(&s1, mpnw)
        mpeq(s0, &c, mpnw)
        var t = mp_real(c[lc...])
        mpeq(s1, &t, mpnw)
        c[lc...] = t[0...]
    } // mpcdiv

    ///   Sets the MPC number B equal to A.
    static func mpceq (_ a:mp_real, _ b: inout mp_real, _ mpnw:Int) {
        let la = a[0]
        let lb = b[0]
        if (mpnw < 4 || a[0] < abs(a[2]) + 4 || a[la] < abs(a[la+2]) + 4 || b[0] < mpnw + 6 || b[lb] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        var ia = sign(1, a[2])
        var na = min(Int(abs(a[2])), mpnw)
        if na == 0  {
//            b[1] = mpnw
//            b[2] = 0
//            b[3] = 0
//            b[4] = 0
//            b[5] = 0
            mp_zero(&b, mpnw)
//            goto 110
        } else {
            b[1] = mpnw
            b[2] = sign(na, ia)
            b[3...na+3] = a[3...na+3]
//            for i in 2...na+2 {
//                b[i+1] = a[i+1]
//            }
//
            b[na+4] = 0
            b[na+5] = 0
        }
        
        // 110 continue
        ia = sign(1, a[la+2])
        na = min(Int(abs(a[la+2])), mpnw)
        if na == 0  {
            b[lb+1] = mpnw
            b[lb+2] = 0
            b[lb+3] = 0
            b[lb+4] = 0
            b[lb+5] = 0
            return
        }
        b[lb+1] = mpnw
        b[lb+2] = sign(na, ia)
        b[lb+3...na+lb+3] = a[la+3...na+la+3]
//        for i in 2...na+2 {
//            b[i+lb+1] = a[i+la+1]
//        }
//
        b[na+lb+4] = 0
        b[na+lb+5] = 0
    } // mpceq

    ///   This routine multiplies the MPC numbers A and B.
    static func mpcmul (_ a:mp_real, _ b: mp_real, _ c: inout mp_real, _ mpnw:Int) {
        var s0 = mp_init(size: mpnw+7), s1 = mp_init(size: mpnw+7), s2 = mp_init(size: mpnw+7)
        var s3 = mp_init(size: mpnw+7)
        let la = a[0]
        let lb = b[0]
        let lc = c[0]
        if (mpnw < 4 || a[0] < abs (a[2]) + 4 || a[la] < abs (a[la+2]) + 4 ||
            b[0] < abs (b[2]) + 4 || b[lb] < abs (b[lb+2]) + 4 ||
            c[0] < mpnw + 6 || c[lc] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        let mpnw1 = mpnw + 1
        s0[0] = mpnw + 7
        s1[0] = mpnw + 7
        s2[0] = mpnw + 7
        s3[0] = mpnw + 7
        
        mpmul(a, b, &s0, mpnw1)
        mpmul(mp_real(a[la...]), mp_real(b[lb...]), &s1, mpnw1)
        mpsub(s0, s1, &s2, mpnw1)
        mpmul(a, mp_real(b[lb...]), &s0, mpnw1)
        mpmul(mp_real(a[la...]), b, &s1, mpnw1)
        mpadd(s0, s1, &s3, mpnw1)
        
        mproun(&s2, mpnw)
        mproun(&s3, mpnw)
        mpeq(s2, &c, mpnw)
        var t = mp_real(c[lc...])
        mpeq(s3, &t, mpnw)
        c[lc...] = t[0...]
    } // mpcmul

    ///   This computes the N-th power of the MPC number A and returns the MPC result
    ///   in B.  When N is zero, 1 is returned.  When N is negative, the reciprocal
    ///   of A ^ |N| is returned.
    static func mpcnpwr (_ a:mp_real, _ n:Int, _ b: inout mp_real, _ mpnw:Int) {
        //    implicit none
        //    integer, intent(in):: mpnw, n
        //    integer j, kk, kn, la, lb, lc, mn, mpnw1, na, nn
        //    real (mprknd) cl2, t1
        var s0 = mp_init(size: 2*mpnw+14), s1 = mp_init(size: 2*mpnw+14), s2 = mp_init(size: 2*mpnw+14)
        let cl2 = 1.4426950408889633
        //    integer (mpiknd), intent(in):: a(0:)
        //    integer (mpiknd), intent(out):: b(0:)
        //    integer (mpiknd) s0(0:2*mpnw+13), s1(0:2*mpnw+13), s2(0:2*mpnw+13)
        
        let la = a[0]
        let lb = b[0]
        if (mpnw < 4 || a[0] < abs(a[2]) + 4 || a[la] < abs(a[la+2]) + 4 || b[0] < mpnw + 6 || b[lb] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        let na = min(max(Int(abs(a[2])), Int(abs(a[la+2]))), mpnw)
        if na == 0 {
            if n >= 0 {
//                b[1] = mpnw
//                b[2] = 0
//                b[3] = 0
//                b[4] = 0
//                b[5] = 0
                mp_zero(&b, mpnw)
                return // goto 120
            } else {
                print("*** \(#function): Argument is zero and N is negative or zero.")
                mpabrt (57)
            }
        }
        
        let mpnw1 = mpnw + 1
        let lc = mpnw + 7
        s0[0] = mpnw + 7
        s0[lc] = mpnw + 7
        s1[0] = mpnw + 7
        s1[lc] = mpnw + 7
        s2[0] = mpnw + 7
        s2[lc] = mpnw + 7
        
        let nn = abs(n)
        if (nn == 0) {
            mpdmc(1.0, 0, &b, mpnw)
            var t = mp_real(b[lb...])
            mpdmc(0.0, 0, &t, mpnw)
            b[lb...] = t[0...]
            return
        } else if (nn == 1) {
            mpceq (a, &s2, mpnw1)
            //  goto 110
        } else if (nn == 2) {
            mpcmul (a, a, &s2, mpnw1)
            //  goto 110
        } else {
            ///   Determine the least integer MN such that 2 ^ MN .GT. NN.
            let t1 = Double(nn)
            let mn = Int(cl2 * Double.log(t1) + 1.0 + mprdfz)
            mpdmc(1.0, 0, &s2, mpnw1)
            var t = mp_real(s2[lc...])
            mpdmc(0.0, 0, &t, mpnw1)
            s2[lc...] = t[0...]
            mpceq(a, &s0, mpnw1)
            var kn = nn
            
            ///   Compute B ^ N using the binary rule for exponentiation.
            for j in 1...Int(mn) {
                let kk = kn / 2
                if kn != 2 * kk {
                    mpcmul(s2, s0, &s1, mpnw1)
                    mpceq(s1, &s2, mpnw1)
                }
                kn = kk
                if j < mn {
                    mpcmul(s0, s0, &s1, mpnw1)
                    mpceq(s1, &s0, mpnw1)
                }
            }
            
            
        }
        // 110 continue
        ///   Compute reciprocal if N is negative.
        if n < 0 {
            mpdmc(1.0, 0, &s1, mpnw1)
            var t = mp_real(s1[lc...])
            mpdmc(0.0, 0, &t, mpnw1)
            s1[lc...] = t[0...]
            mpcdiv(s1, s2, &s0, mpnw1)
            mpceq(s0, &s2, mpnw1)
        }
        
        ///   Restore original precision level.
        mproun(&s2, mpnw)
        var t = mp_real(s2[lc...])
        mproun(&t, mpnw)
        s2[lc...] = t[0...]
        mpceq(s2, &b, mpnw)
    } // mpcnpwr

    ///   This routine returns the conjugate of the MPC argument A.
    static func mpconjg(_ a:mp_real, _ b: inout mp_real, _ mpnw:Int) {
        let la = a[0]
        let lb = b[0]
        if (mpnw < 4 || a[0] < abs(a[2]) + 4 || a[la] < abs(a[la+2]) + 4 || b[0] < mpnw + 6 || b[lb] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        mpceq(a, &b, mpnw)
        b[lb+2] = -b[lb+2]
    } // mpconjg

    ///   This routine returns the square root of the MPC argument A.
    ///   The formula is:
    ///
    ///   1/Sqrt[2] * (Sqrt[r + a1] + I * a2 / Sqrt[r + a1])  if a1 >= 0, or
    ///   1/Sqrt[2] * (|a2| / Sqrt[r - a1] + I * Sgn[a2] * Sqrt[r - a1]) if a1 < 0,
    ///
    ///   where r = Sqrt[a1^2 + a2^2], and a1 and a2 are the real and imaginary
    ///   parts of A.
    static func mpcsqrt (_ a:mp_real, _ b: inout mp_real, _ mpnw:Int) {
        var s0 = mp_init(size: mpnw+7), s1 = mp_init(size: mpnw+7), s2 = mp_init(size: mpnw+7)
        var s3 = mp_init(size: mpnw+7), s4 = mp_init(size: mpnw+7)
        let la = a[0]
        let lb = b[0]
        if (mpnw < 4 || a[0] < abs (a[2]) + 4 || a[la] < abs (a[la+2]) + 4 || b[0] < mpnw + 6 || b[lb] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        let mpnw1 = mpnw + 1
        s0[0] = mpnw + 7
        s1[0] = mpnw + 7
        s2[0] = mpnw + 7
        s3[0] = mpnw + 7
        s4[0] = mpnw + 7
        
        mpmul (a, a, &s0, mpnw1)
        mpmul (mp_real(a[la...]), mp_real(a[la...]), &s1, mpnw1)
        mpadd (s0, s1, &s2, mpnw1)
        mpsqrt (s2, &s0, mpnw1)
        
        if (a[2] >= 0) {
            mpadd (s0, a, &s1, mpnw1)
            mpsqrt (s1, &s0, mpnw1)
            mpdiv (mp_real(a[la...]), s0, &s1, mpnw1)
        } else {
            mpsub (s0, a, &s2, mpnw1)
            mpsqrt (s2, &s1, mpnw1)
            mpdiv (mp_real(a[la...]), s1, &s0, mpnw1)
            s0[2] = abs(s0[2])
            s1[2] = sign(1, a[la+2])
        }
        
        mpdmc (0.5, 0, &s3, mpnw1)
        mpsqrt (s3, &s2, mpnw1)
        mpmul (s0, s2, &s3, mpnw1)
        mpmul (s1, s2, &s4, mpnw1)
        
        mproun (&s3, mpnw)
        mproun (&s4, mpnw)
        mpeq (s3, &b, mpnw)
        var t = mp_real(b[lb...])
        mpeq (s4, &t, mpnw)
        b[lb...] = t[0...]
    } // mpcsqrt

    ///   This routine subtracts the MPC numbers A and B.
    static func mpcsub (_ a:mp_real, _ b: mp_real, _ c: inout mp_real, _ mpnw:Int) {
        let la = a[0]
        let lb = b[0]
        let lc = c[0]
        if (mpnw < 4 || a[0] < abs(a[2]) + 4 || a[la] < abs(a[la+2]) + 4
            || b[0] < abs(b[2]) + 4 || b[lb] < abs(b[lb+2]) + 4 ||
            c[0] < mpnw + 6 || c[lc] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        mpsub (a, b, &c, mpnw)
        var t = mp_real(c[lc...])
        mpsub (mp_real(a[la...]), mp_real(b[lb...]), &t, mpnw)
        c[lc...] = t[0...]
    } // mpcsub

    // MARK: - More elementary operations
    ///   This routine compares the MPR numbers A and B and returns in IC the value
    ///   -1, 0, or 1 depending on whether A < B, A = B, or A > B.
    ///   Note that the first and second words do NOT need to be the same for the
    ///   result to be "equal".
    static func mpcpr (_ a:mp_real, _ b: mp_real, _ ic: inout Int, _ mpnw:Int) {
        var s0 = mp_init(size: mpnw+6)
        
        if (mpnw < 4 || a[0] < abs(a[2]) + 4 || b[0] < abs(b[2]) + 4) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        s0[0] = mpnw + 6
        mpsub (a, b, &s0, mpnw)
        if (s0[2] < 0) {
            ic = -1
        } else if (s0[2] == 0) {
            ic = 0
        } else {
            ic = 1
        }
    } // mpcpr

    ///   This divides A by B and returns the result in C.
    ///   This function employs the following Newton-Raphson iteration, which
    ///   converges to 1 / B:
    ///
    ///    X_{k+1} = X_k + (1 - X_k \* B) \* X_k
    ///
    ///   where the multiplication () \* X_k is performed with only half of the
    ///   normal level of precision.  These iterations are performed with a
    ///   working precision level MPNW that is dynamically changed, approximately
    ///   doubling with each iteration (except that at iteration NIT before the
    ///   final iteration, the iteration is repeated without doubling the
    ///   precision, in order to enhance accuracy).  The final iteration is
    ///   performed as follows (this is due to A. Karp):
    ///
    ///    A / B = (A \* X_n) + [A - (A \* X_n) \* B] \* X_n  (approx.)
    ///
    ///   where the multiplications A \* X_n and [] \* X_n are performed with only
    ///   half of the final level of precision.
    static func mpdiv (_ a:mp_real, _ b: mp_real, _ c: inout mp_real, _ mpnw:Int) {
        let cl2 = 1.4426950408889633
        let nit = 3
        var s0 = mp_init(size: mpnw+7), s1 = mp_init(size: mpnw+7), s2 = mp_init(size: mpnw+7)
        var s3 = mp_init(size: mpnw+7)
        if (mpnw < 4 || a[0] < abs(a[2]) + 4 || b[0] < abs(b[2]) + 4 || c[0] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        let na = min(Int(abs(a[2])), mpnw)
        let nb = min(Int(abs(b[2])), mpnw)
        
        if (na == 0) {
//            c[1] = mpnw
//            c[2] = 0
//            c[3] = 0
//            c[4] = 0
//            c[5] = 0
            mp_zero(&c, mpnw)
            return
        }
        if nb == 0 {
            print("*** \(#function): Divisor is zero.")
            mpabrt (33)
            return
        }
        
        s0[0] = mpnw + 7
        s1[0] = mpnw + 7
        s2[0] = mpnw + 7
        s3[0] = mpnw + 7
        
        ///   Determine the least integer MQ such that 2 ^ MQ .GE. MPNW.
        var t1 = Double(mpnw)
        let mq = Int(cl2 * Double.log(t1) + 1.0 - mprdfz)
        
        ///   Compute the initial approximation of 1 / B.
        var n = 0
        mpmdc (b, &t1, &n, mpnw)
        let t2 = 1.0 / t1
        mpdmc(t2, -n, &s2, mpnw)
        mpdmc(1.0, 0, &s3, mpnw)
        var mpnw1 = 5
        var iq = 0
        var nw1 = mpnw1
        var nw2 = mpnw1
        
        ///   Perform the Newton-Raphson iteration described above with a dynamically
        ///   changing precision level MPNW (one greater than powers of two).
        for k in 1..<mq {
            if (k > 2) {
                nw1 = mpnw1
                mpnw1 = min (2 * mpnw1 - 2, mpnw) + 1
                nw2 = mpnw1
            }
            
            // 100  continue
            while true {
                mpmul (b, s2, &s1, nw2)
                mpsub (s3, s1, &s0, nw2)
                mpmul (s2, s0, &s1, nw1)
                mpadd (s2, s1, &s0, nw2)
                mpeq (s0, &s2, nw2)
                if k == mq - nit && iq == 0 { iq = 1; continue }
                break
            }
        }
        
        ///   Perform last iteration using Karp"s trick.
        nw1 = mpnw1
        mpnw1 = min (2 * mpnw1 - 1, mpnw) + 1
        nw2 = mpnw1
        
        mpmul (a, s2, &s0, nw1)
        mpmul (s0, b, &s1, nw2)
        mpsub (a, s1, &s3, nw2)
        mpmul (s3, s2, &s1, nw1)
        mpadd (s0, s1, &s2, nw2)
        
        ///   Restore original precision level.
        mproun (&s2, mpnw)
        mpeq (s2, &c, mpnw)
    } // mpdiv
    
    ///   This routine divides the MPR number A by the DP number B to yield C.
    ///   NOTE however that the product is not fully accurate unless B is an exact
    ///   binary value.
    ///   Examples of exact binary values (good): 123456789.0, 0.25d0, -5.3125d0.
    ///   Examples of inexact binary values (bad): 0.1d0, 123467.8d0, -3333.3d0.
    static func mpdivd (_ a:mp_real, _ b:Double, _ c: inout mp_real, _ mpnw: Int) {
        //    implicit none
        //    integer, intent(in):: mpnw
        //    integer i, ia, ib, j, k, na, nbth, n1
        let nbth = mpnbt / 2
        //    real (mprknd), intent(in):: b
        //    real (mprknd) bb, bdh, bdvd, rdh
        let bdh = Double.pow(2.0, Double(nbth))
        let rdh = Double.pow(0.5, Double(nbth))
        var cc = mp_init(size: 2*mpnw+11), d = mp_init(size: 2*mpnw+11)
        //    integer (mpiknd), intent(in):: a(0:)
        //    integer (mpiknd), intent(out):: c(0:)
        //    integer (mpiknd) cc(0:2*mpnw+10), d(0:2*mpnw+10), ibb, &
        //    b1, b2, c11, c12, c21, c22, d1, d2, td, t1, t2, t3
        if (mpnw < 4 || a[0] < abs (a[2]) + 4 || c[0] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        ///   Check for zero inputs.
        let ia = sign(1, a[2])
        let na = min(Int(abs(a[2])), mpnw)
        let ib = sign(1, Int(b))
        if na == 0 {
//            c[1] = mpnw
//            c[2] = 0
//            c[3] = 0
//            c[4] = 0
//            c[5] = 0
            mp_zero(&c, mpnw)
            return
        } else if b == 0.0 {
            print("*** \(#function): Divisor is zero.")
            mpabrt (31)
        } else if b == 1.0 {
            mpeq(a, &c, mpnw)
            return
        }
        var bb = abs(b)
        var n1 = 0
        
        ///   Reduce BB to within 1 and MPBDX.
        if bb >= Double(mpbdx) {
            for k in 1...100 {
                bb = bb / Double(mpbdx)
                if bb < Double(mpbdx) {
                    n1 = n1 + k
                    break // goto 120
                }
            }
        } else if bb < 1.0 {
            for k in 1...100 {
                bb = Double(mpbdx) * bb
                if bb >= 1.0 {
                    n1 = n1 - k
                    break // goto 120
                }
            }
        }
        
        // 120  continue
        let ibb = Int(bb)
        
        ///   If bb is not an integer, mpdiv instead.
        if (bb != Double(ibb)) {
            d[0] = mpnw + 6
            d[1] = mpnw
            mpdmc (b, 0, &d, mpnw)
            mpdiv (a, d, &c, mpnw)
            return
        }
        
        cc[0] = 0
        cc[1] = mpnw
        cc[2] = sign(mpnw, ia * ib)
        cc[3...2*mpnw+10] = [Int](repeating: 0, count: 2*mpnw+8)[0...]
        
        ///   Split D array into half-word chunks.
        d[0...3] = [0, 0, 0, 0]
        for i in 0...na {
            let c11 = a[i+4] >> nbth
            let c12 = a[i+4] - (c11 << nbth)
            d[2*i+4] = c11
            d[2*i+5] = c12
        }
        
        d[2*na+6...2*mpnw+10] = [Int](repeating: 0, count: 2*mpnw+5-2*na)[0...]
        let b1 = ibb >> nbth
        let b2 = ibb - (b1 << nbth)
        
        ///   Perform short division algorithm, after splitting inputs.
        ///   Double precision is employed to find and refine the trial divisor.
        var t1:Int
        for j in 3...2*mpnw+5 {
            let bdvd = bdh * Double(d[j]) + Double(d[j+1]) + Double(d[j+2]) * rdh
            let td = Int(floor(bdvd / bb))
            t1 = b1 * td
            let c11 = t1 >> nbth
            let c12 = t1 - (c11 << nbth)
            let t2 = b2 * td
            let c21 = t2 >> nbth
            let c22 = t2 - (c21 << nbth)
            let d1 = c12 + c21 + (c11 << nbth)
            let d2 = c22
            d[j] = d[j] - d1
            d[j+1] = d[j+1] - d2 + (d[j] << nbth)
            cc[j+1] = td
        }
        
        ///  Release carries on the full cc vector.
        t1 = 0
        for i in stride(from: 2*mpnw+5, through: 3, by: -1) {
            let t3 = t1 + cc[i+1]
            t1 = t3 >> nbth
            cc[i+1] = t3 - (t1 << nbth)
        }
        
        cc[3] = cc[3] + t1
        
        ///  Recombine half words into full words.
        c[1] = mpnw
        c[2] = sign(mpnw, ia * ib)
        c[3] = cc[3]
        c[4] = cc[4]
        
        for i in 0...mpnw+1 {
            c[i+4] = cc[2*i+4] << nbth + cc[2*i+5]
        }
        
        ///   If c[3] is nonzero, shift the result one cell right.
        if c[3] != 0 {
            n1 = n1 + 1
            c[2] = sign(abs(c[2]) + 1, c[2])
            
            for i in stride(from: mpnw+4, through: 3, by: -1) {
                c[i+1] = c[i]
            }
        }
        c[3] = a[3] - n1
        mproun(&c, mpnw)
    } // mpdivd
    
    static func mpmask13(_ b:Double) -> Double {
        let b13x = Double.pow(2.0, 13.0)
        let t1 = b13x * abs(b)
        return abs(abs(b) + t1) - abs(t1)
    }
    
    static func mpmask23(_ b:Double) -> Double {
        let b23x = Double.pow(2.0, 23.0)
        let t1 = b23x * abs(b)
        return abs(abs(b) + t1) - abs(t1)
    }

    ///   This routine divides the MPR number A by the DP number B to yield C.
    ///   In contrast to mpdivd, this routine only allows 40 significant bits
    ///   (approximately 12 significant decimal digits) in B.  If more nonzero bits
    ///   are present in B (likely due to inexact binary value), an error is flagged.
    ///   Examples of exact binary values (good): 123456789.0, 0.25d0, -5.3125d0.
    ///   Examples of inexact binary values (bad): 0.1d0, 123467.8d0, -3333.3d0.
    static func mpdivd40 (_ a:mp_real, _ b:Double, _ c: inout mp_real, _ mpnw:Int) {
        if (mpnw < 4 || a[0] < abs(a[2]) + 4 || c[0] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        ///   Check whether B has more than 40 significant bits (actually whether the
        ///   trailing 13 bits are zero).
        let t2 = mpmask13(b)
        if t2 == abs(b) {
            mpdivd (a, b, &c, mpnw)
        } else {
            print("*** \(#function): DP value has more than 40 significant bits:",
                  "\(b) and thus very likely represents an unintended loss of accuracy.",
                  "Fix the issue, or else use functions mpprod, mpquot, mpreald or mpcmplxdc.",
                  "See documentation for details.")
            mpabrt (81)
        }
    } // mpdivd40

    // MARK: - To/From MPReal number conversions
    
    ///   This routine converts the DP number A * 2^N to MPR form in B.
    ///   NOTE however that the product is not fully accurate unless A is an exact
    ///   binary value.
    ///   Examples of exact binary values (good): 123456789.0, 0.25d0, -5.3125d0.
    ///   Examples of inexact binary values (bad): 0.1d0, 123467.8d0, -3333.3d0.
    static func mpdmc (_ a:Double, _ n:Int, _ b: inout mp_real, _ mpnw:Int) {
        if (mpnw < 4 || b[0] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        ///   Check for zero.
        if a == 0.0 {
//            b[1] = mpnw
//            b[2] = 0
//            b[3] = 0
//            b[4] = 0
//            b[5] = 0
            mp_zero(&b, mpnw)
            return
        }
        var n1 = n / mpnbt
        let n2 = n - mpnbt * n1
        var aa = abs(a) * Double.pow(2.0, Double(n2))
        
        ///   Reduce AA to within 1 and MPBDX.
        if aa >= Double(mpbdx) {
            for k in 1...100 {
                aa = aa / Double(mpbdx)
                if aa < Double(mpbdx) {
                    n1 = n1 + k
                    break // goto 120
                }
            }
            
        } else if aa < 1.0 {
            for k in 1...100 {
                aa = aa * Double(mpbdx)
                if aa >= 1.0 {
                    n1 = n1 - k
                    break // goto 120
                }
            }
        }
        
        ///   Store successive sections of AA into B.
        b[3] = n1
        b[4] = Int(aa)
        aa = Double(mpbdx) * (aa - Double(b[4]))
        b[5] = Int(aa)
        b[6] = 0
        b[7] = 0
        b[8] = 0
        
        var ii = 3
        for i in stride(from: 6, through: 3, by: -1) {
            if b[i+1] != 0 { ii = i; break /*goto 140*/ }
        }
        b[1] = mpnw
        aa = Double(ii - 2)
        b[2] = Int(sign(aa, a))
    } // mpdmc

    ///   This routine converts the DP number A \* 2^N to MPR form in B.
    ///   In contrast to mpdmc, this routine only allows 40 significant bits
    ///   (approximately 12 significant decimal digits) in A.  If more nonzero bits
    ///   are present in A (likely due to inexact binary value), an error is flagged.
    ///   Examples of exact binary values (good): 123456789.0, 0.25, -5.3125.
    ///   Examples of inexact binary values (bad): 0.1, 123467.8, -3333.3.
    static func mpdmc40 (_ a:Double, _ n:Int, _ b: inout mp_real, _ mpnw:Int) {
        if (mpnw < 4 || b[0] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        ///   Check whether A has more than 40 significant bits (actually whether
        ///   the trailing 13 bits are zero).
        let t2 = mpmask13(a)
        if t2 == abs(a) {
            mpdmc(a, n, &b, mpnw)
        } else {
            print("*** \(#function): DP value has more than 40 significant bits:",
                  "\(a) and thus very likely represents an unintended loss of accuracy.",
                  "Fix the issue, or else use functions mpprod, mpquot, mpreald or mprealdm.",
                  "See documentation for details.")
            mpabrt (82)
        }
    } // mpdmc40

    ///   Sets the MPR number B equal to the MPR number A.
    static func mpeq (_ a: mp_real, _ b: inout mp_real, _ mpnw:Int) {
        if (mpnw < 4 || a[0] < abs (a[2]) + 4 || b[0] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        let ia = sign(1, a[2])
        let na = min(Int(abs(a[2])), mpnw)
        if na == 0  {
//            b[1] = mpnw
//            b[2] = 0
//            b[3] = 0
//            b[4] = 0
//            b[5] = 0
            mp_zero(&b, mpnw)
            return
        }
        b[1] = mpnw
        b[2] = sign (na, ia)
        b[3...na+3] = a[3...na+3]
        //          for i = 2, na + 2
        //            b(i+1) = a(i+1)
        //          }
        b[na+4] = 0
        b[na+5] = 0
        // 110 continue
    } // mpeq

    ///   Sets B to the integer part of the MPR number A and sets C equal to the
    ///   fractional part of A.  Note this is NOT the quite same as the greatest
    ///   integer function as often defined in some mathematical books and papers.
    ///   Examples:  If A = 1.75, then B = 1.0, C = 0.75.
    ///     If A = -3.25, then B = -3.0, C = -0.25.
    static func mpinfr (_ a:mp_real, _ b: inout mp_real, _ c: inout mp_real, _ mpnw:Int) {
        if (mpnw < 4 || a[0] < abs (a[2]) + 4 || b[0] < mpnw + 6 || c[0] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        ///   Check if  A  is zero.
        let ia = sign(1, a[2])
        let na = min(Int(abs(a[2])), mpnw)
        let ma = a[3]
        if na == 0  {
//            b[1] = mpnw
//            b[2] = 0
//            b[3] = 0
//            b[4] = 0
//            b[5] = 0
            mp_zero(&b, mpnw)
//            c[1] = mpnw
//            c[2] = 0
//            c[3] = 0
//            c[4] = 0
//            c[5] = 0
            mp_zero(&c, mpnw)
            return // goto 120
        }
        
        if (ma >= mpnw - 1) {
            print("*** \(#function): Argument is too large.")
            mpabrt (40)
        }
        
        ///   Place integer part in  B.
        let nb = min(max(ma + 1, 0), na)
        if nb == 0 {
//            b[1] = mpnw
//            b[2] = 0
//            b[3] = 0
//            b[4] = 0
//            b[5] = 0
            mp_zero(&b, mpnw)
        } else {
            b[1] = mpnw
            b[2] = sign(nb, ia)
            b[3] = ma
            b[nb+4] = 0
            b[nb+5] = 0
            b[4...nb+3] = a[4...nb+3]
//            for i = 3, nb + 2 {
//                b(i+1) = a(i+1)
//            }
        }
        
        ///   Place fractional part in C.
        let nc = na - nb
        if nc <= 0 {
//            c[1] = mpnw
//            c[2] = 0
//            c[3] = 0
//            c[4] = 0
//            c[5] = 0
            mp_zero(&c, mpnw)
        } else {
            c[1] = mpnw
            c[2] = sign(nc, ia)
            c[3] = ma - nb
            c[nc+4] = 0
            c[nc+5] = 0
            c[4...nc+3] = a[nb+4...nc+nb+3]
//            for i = 3, nc + 2 {
//                c(i+1) = a(i+nb+1)
//            }
        }
        
        ///   Fix up results.  B may have trailing zeros and C may have leading zeros.
        mproun(&b, mpnw)
        mproun(&c, mpnw)
    } // mpinfr

    ///   This returns a DP approximation the MPR number A in the form B * 2^n.
    static func mpmdc (_ a:mp_real, _ b: inout Double, _ n: inout Int, _ mpnw:Int) {
        if (mpnw < 4 || a[0] < abs (a[2]) + 4) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        if a[2] == 0  {
            b = 0.0
            n = 0
            return
        }
        
        var na = abs(a[2])
        var aa = Double(a[4])
        if na >= 2 { aa = aa + Double(a[5]) / Double(mpbdx) }
        n = mpnbt * a[3]
        b = sign(aa, Double(a[2]))
        
        ///   Reduce b to within 1 and 2.
        na = Int(Double.log(abs(b)) / Double.log(2.0) + Double(mprdfz))
        b = b / Double.pow(2.0, Double(na))
        n = n + na
        if abs(b) < 1.0 {
            b = 2.0 * b
            n = n - 1
        } else if abs(b) > 2.0 {
            b = 0.5 * b
            n = n + 1
        }
    } // mpmdc

    ///   This routine multiplies MPR numbers A and B to yield C.
    
    ///   This routine returns up to MPNW mantissa words of the product.  If the
    ///   complete double-long product of A and B is desired (for example in large
    ///   integer applications), then MPNW must be at least as large as the sum of
    ///   the mantissa lengths of A and B.  In other words, if the precision levels
    ///   of A and B are both 64 words, then MPNW must be at least 128 words to
    ///   produce the complete double-long product in C.
    static func mpmul(_ a:mp_real, _ b: mp_real, _ c: inout mp_real, _ mpnw:Int) {
        let nbth = mpnbt / 2
        var d = mp_init(size: mpnw+7)
        
        if (mpnw < 4 || a[0] < abs (a[2]) + 4 || b[0] < abs (b[2]) + 4 || c[0] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        let ia = sign(1, a[2])
        let ib = sign(1, b[2])
        let na = min(Int(abs (a[2])), mpnw)
        let nb = min(Int(abs (b[2])), mpnw)
        var nc = min(na + nb, mpnw)
        
        if na == 0 || nb == 0 {
            ///   One of the inputs is zero -- result is zero.
//            c[1] = mpnw
//            c[2] = 0
//            c[3] = 0
//            c[4] = 0
//            c[5] = 0
            mp_zero(&c, mpnw)
            return
        }
        
        if na == 1 && a[4] == 1 {
            ///   A is 1 or -1 -- result is B or -B.
            c[1] = mpnw
            c[2] = sign(nb, ia * ib)
            c[3] = a[3] + b[3]
            c[4...nb+3] = b[4...nb+3]
//            for i = 3, nb + 2 {
//                c(i+1) = b(i+1)
//            }
            c[nb+4] = 0
            c[nb+5] = 0
            return
        } else if nb == 1 && b[4] == 1 {
            ///   B is 1 or -1 -- result is A or -A.
            c[1] = mpnw
            c[2] = sign (na, ia * ib)
            c[3] = a[3] + b[3]
            c[4...na+3] = a[4...na+3]
//            for i = 3, na + 2 {
//                c(i+1) = a(i+1)
//            }
            c[na+4] = 0
            c[na+5] = 0
            return
        }
        
        ///   For very high precision, mpmulx.
        if (na > mpmlxm && nb > mpmlxm) {
            mpmulx(a, b, &c, mpnw)
            return
        }
        
        var dd = a[3] + b[3]
        d[0] = mpnw + 6
        d[1] = mpnw
        d[2] = sign(nc, ia * ib)
        d[3...nc+6] = [Int](repeating: 0, count: nc+4)[0...]
//        for i = 2, nc + 5 {
//            d(i+1) = 0
//        }
        
        ///   Perform ordinary long multiplication algorithm, after splitting inputs.
        ///   Accumulate at most MPNW+2 mantissa words of the product.
        var t1:Int
        for j in 3...na+2 {
            let j3 = j - 3
            let n2 = min(nb + 2, mpnw + 4 - j3)
            let a1 = a[j+1] >> nbth
            let a2 = a[j+1] - (a1 << nbth)
            
            for i in 3...n2 {
                let b1 = b[i+1] >> nbth
                let b2 = b[i+1] - (b1 << nbth)
                let c1 = a1 * b2 + a2 * b1
                let c2 = c1 >> nbth
                let c3 = c1 - (c2 << nbth)
                d[i+j3] = d[i+j3] + a1 * b1 + c2
                d[i+j3+1] = d[i+j3+1] + a2 * b2 + (c3 << nbth)
            }
            
            ///  Release carries on the just-computed section of the d vector.
            t1 = 0
            for i in stride(from: n2, through: 3, by: -1) {
                let t3 = t1 + d[i+j3+1]
                t1 = t3 >> mpnbt
                d[i+j3+1] = t3 - (t1 << mpnbt)
            }
            d[j3+3] = d[j3+3] + t1
        }
        
        ///  Release carries on the full d vector.
        t1 = 0
        for i in stride(from: nc+1, through: 1, by: -1) {
            let t3 = t1 + d[i+3]
            t1 = t3 >> mpnbt
            d[i+3] = t3 - (t1 << mpnbt)
        }
        d[3] = d[3] + t1
        
        ///   If d[3] is nonzero, shift the result one cell right.
        if (d[3] != 0) {
            dd = dd + 1
            nc = min (nc + 1, mpnw)
            d[2] = sign (nc, Int(d[2]))
            
            for i in stride(from: nc+4, through: 3, by: -1) {
                d[i+1] = d[i]
            }
        }
        
        d[3] = dd
        c[1...nc+5] = d[1...nc+5]
//        for i in 1, nc + 5 {
//            c(i) = d(i)
//        }
        mproun(&c, mpnw)
    } // mpmul

    ///   This routine multiplies the MPR number A by the DP number B to yield C.
    ///   NOTE however that the product is not fully accurate unless B is an exact
    ///   binary value.
    ///   Examples of exact binary values (good): 123456789.0, 0.25, -5.3125.
    ///   Examples of inexact binary values (bad): 0.1, 123467.8, -3333.3.
    static func mpmuld (_ a:mp_real, _ b:Double, _ c: inout mp_real, _ mpnw:Int) {
        let nbth = mpnbt / 2
        var d = mp_init(size: mpnw+7)
        
        if (mpnw < 4 || a[0] < abs(a[2]) + 4 || c[0] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        ///   Check for zero inputs.
        let ia = sign(1, a[2])
        let na = min(Int(abs(a[2])), mpnw)
        let ib = sign(1.0, b)
        if na == 0 || b == 0.0 {
//            c[1] = mpnw
//            c[2] = 0
//            c[3] = 0
//            c[4] = 0
//            c[5] = 0
            mp_zero(&c, mpnw)
            return
        } else if b == 1.0 {
            mpeq(a, &c, mpnw)
            return
        }
        
        var bb = abs(b)
        var n1 = 0
        
        ///   Reduce BB to within 1 and MPBDX.
        if bb >= Double(mpbdx) {
            for k in 1...100 {
                bb = bb / Double(mpbdx)
                if bb < Double(mpbdx) {
                    n1 = n1 + k
                    break // goto 120
                }
            }
        } else if bb < 1.0 {
            for k in 1...100 {
                bb = Double(mpbdx) * bb
                if bb >= 1.0 {
                    n1 = n1 - k
                    break // goto 120
                }
            }
        }
        
       //  120  continue
        
        let ibb = Int(bb)
        
        ///   If bb is not an integer, mpmul instead.
        if bb != Double(ibb) {
            d[0] = mpnw + 6
            d[1] = mpnw
            mpdmc(b, 0, &d, mpnw)
            mpmul(a, d, &c, mpnw)
            return
        }
        
        d[0] = mpnw + 6
        d[1] = mpnw
        d[2] = sign(na, ia * Int(ib))
        d[3...na+6] = [Int](repeating: 0, count: na+4)[0...]
//        for i in 2...na+5 {
//            d[i+1] = 0
//        }
        
        let b1 = ibb >> nbth
        let b2 = ibb - (b1 << nbth)
        
        ///   Perform short multiplication algorithm, after splitting inputs.
        for j in 3...na+3 {
            let a1 = a[j+1] >> nbth
            let a2 = a[j+1] - (a1 << nbth)
            let c1 = a1 * b2 + a2 * b1
            let c2 = c1 >> nbth
            let c3 = c1 - (c2 << nbth)
            d[j] = d[j] + a1 * b1 + c2
            d[j+1] = d[j+1] + a2 * b2 + (c3 << nbth)
        }
        
        ///  Release carries on the full d vector.
        var t1 = 0
        for i in stride(from: na+3, through: 3, by: -1) {
            let t3 = t1 + d[i+1]
            t1 = t3 >> mpnbt
            d[i+1] = t3 - (t1 << mpnbt)
        }
        
        d[3] = d[3] + t1
        
        ///   If d[3] is nonzero, shift the result one cell right.
        if d[3] != 0 {
            n1 = n1 + 1
            d[2] = sign(abs(d[2]) + 1, d[2])
            
            for i in stride(from: na+4, through: 3, by: -1) {
                d[i+1] = d[i]
            }
        }
        
        d[3] = a[3] + n1
        
        ///   Copy d to c and round.
        c[1...na+5] = d[1...na+5]
//        for i = 1, na + 5 {
//            c(i) = d(i)
//        }
//
        mproun (&c, mpnw)
    } // mpmuld
    
    ///   This routine multiples the MP number A by the DP number B to yield C.
    ///   In contrast to mpmuld, this routine only allows 40 significant bits
    ///   (approximately 12 significant decimal digits) in B.  If more nonzero bits
    ///   are present in B (likely due to inexact binary value), an error is flagged.
    ///   Examples of exact binary values (good): 123456789.0, 0.25, -5.3125.
    ///   Examples of inexact binary values (bad): 0.1, 123467.8, -3333.3.
    static func mpmuld40 (_ a:mp_real, _ b:Double, _ c: inout mp_real, _ mpnw:Int) {
        if (mpnw < 4 || a[0] < abs (a[2]) + 4 || c[0] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        ///   Check whether B has more than 40 significant bits (actually whether
        ///   the trailing 13 bits are zero).
        let t2 = mpmask13(b)
        if t2 == abs(b) {
            mpmuld (a, b, &c, mpnw)
        } else {
            print("*** \(#function): DP value has more than 40 significant bits: \(b)",
                      " and thus very likely represents an unintended loss of accuracy.",
                      "Fix the issue, or else use functions mpprod, mpquot, mpreald or mpcmplxdc.",
                      "See documentation for details.")
            mpabrt(83)
        }
    } // mpmuld40

    ///   This routine sets rb = negation of ra.
    static func mpneg ( _ ra:mp_real, _ rb: inout mp_real, _ mpnw:Int) {
        mpeq(ra, &rb, mpnw)
        let na = min(abs(Int(ra[2])), mpnw)
        rb[2] = -sign(na, Int(ra[2]))
    } // mpneg

    ///   This sets B to the nearest integer to the MPR number A.
    ///   Examples:  If A = 1.49, B = 1.0; if A = 3.5, B = 4.0; if A = -2.5, B = -3.0
    static func mpnint( _ a:mp_real, _ b: inout mp_real, _ mpnw:Int) {
        var s0 = mp_init(size: mpnw+7), s1 = mp_init(size: mpnw+7)
        
        if (mpnw < 4 || a[0] < abs (a[2]) + 4 || b[0] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        let ia = sign(1, a[2])
        let na = min(Int(abs(a[2])), mpnw)
        let ma = a[3]
        if na == 0 {
            ///   A is zero -- result is zero.
//            b[1] = mpnw
//            b[2] = 0
//            b[3] = 0
//            b[4] = 0
//            b[5] = 0
            mp_zero(&b, mpnw)
            return
        }
        
        if ma >= mpnw {
            ///   A cannot be represented exactly as an integer.
            print("*** \(#function): Argument is too large.")
            mpabrt(56)
        }
        
        ///   Add or subtract 1/2 from the input, depending on its sign, {
        ///   return the greatest integer.
        s0[0] = mpnw + 6
        s1[0] = mpnw + 6
        
        mpdmc(0.5, 0, &s0, mpnw)
        if ia == 1 {
            mpadd (a, s0, &s1, mpnw)
        } else {
            mpsub (a, s0, &s1, mpnw)
        }
        mpinfr (s1, &b, &s0, mpnw)
    } // mpnint

    ///   This converts the MP number in array D to the standard normalized form
    ///   in A.
    ///   MPNORM assumes that two extra mantissa words are input at the end of D.
    ///   This reduces precision loss when it is necessary to shift the result to
    ///   the left. All words up to index A[2] + 5 in A *must* have data, even if 0.
    static func mpnorm (_ d: inout mp_real, _ a: inout mp_real, _ mpnw:Int) {
        if (mpnw < 4 || d[0] < abs (d[2]) + 4 || a[0] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        var ia = sign(1, d[2])
        var na = min(Int(abs(d[2])), mpnw)
        if na == 0  {
//            a[1] = mpnw
//            a[2] = 0
//            a[3] = 0
//            a[4] = 0
//            a[5] = 0
            mp_zero(&a, mpnw)
            return
        }
        let n4 = na + 4
        var a2 = d[3]
        d[3] = 0
        
        // 110 continue
        var gotoActive:Bool // needed to eliminate goto
        repeat {
            gotoActive = false
            var t1 = 0
            for i in stride(from: n4, through: 3, by: -1) {
                let t3 = t1 + d[i+1]
                t1 = t3 >> mpnbt
                d[i+1] = t3 - (t1 << mpnbt)
            }
            
            d[3] = d[3] + t1
            
            if d[3] < 0 {
                ///   D[3] is negative -- negate all words and re-normalize.
                ia = -ia
                d[4] = d[4] + mpbdx * d[3]
                d[3] = 0
                for i in 2...n4 {
                    d[i+1] = -d[i+1]
                }
                
                gotoActive = true // goto 110
            } else if d[3] > 0 {
                ///   The fixup loops above "spilled" a nonzero number into D[3].  Shift the
                ///   entire number right one cell.  The exponent and length of the result
                ///   are increased by one.
                for i in stride(from: n4, through: 3, by: -1) {
                    a[i+1] = d[i]
                }
                na = min(na + 1, mpnw)
                a2 = a2 + 1
            } else {
                for i in 3...n4 {
                    a[i+1] = d[i+1]
                }
            }
        } while gotoActive
        
        ///   Perform rounding and truncation.
        a[1] = mpnw
        a[2] = sign(na, ia)
        a[3] = a2
        mproun (&a, mpnw)
    } // mpnorm

    // MARK: - Power and root operations
    
    ///   This computes the N-th power of the MPR number A and returns the result
    ///   in B.  When N is zero, 1 is returned.  When N is negative, the reciprocal
    ///   of A ^ |N| is returned.
    static func mpnpwr (_ a:mp_real, _ n:Int, _ b: inout mp_real, _ mpnw:Int) {
        let cl2 = 1.4426950408889633
        var s0 = mp_init(size: mpnw+7), s1 = mp_init(size: mpnw+7), s2 = mp_init(size: mpnw+7)
        
        if (mpnw < 4 || a[0] < abs (a[2]) + 4 || b[0] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        let na = min(Int(abs(a[2])), mpnw)
        if na == 0 {
            if n >= 0 {
//                b[1] = mpnw
//                b[2] = 0
//                b[3] = 0
//                b[4] = 0
//                b[5] = 0
                mp_zero(&b, mpnw)
                return
            } else {
                print("*** \(#function): Argument is zero and N is negative or zero.")
                mpabrt (57)
            }
        }
        
        let mpnw1 = mpnw + 1
        s0[0] = mpnw + 7
        s1[0] = mpnw + 7
        s2[0] = mpnw + 7
        
        let nn = abs(n)
        if nn == 0 {
            mpdmc (1.0, 0, &b, mpnw)
            return
        } else if nn == 1 {
            mpeq(a, &s2, mpnw1)
            // goto 110
        } else if nn == 2 {
            mpmul(a, a, &s2, mpnw1)
            // goto 110
        } else {
            ///   Determine the least integer MN such that 2 ^ MN .GT. NN.
            let t1 = Double(nn)
            let mn = Int(cl2 * Double.log(t1) + 1.0 + Double(mprdfz))
            mpdmc(1.0, 0, &s2, mpnw1)
            mpeq(a, &s0, mpnw1)
            var kn = nn
            
            ///   Compute B ^ N using the binary rule for exponentiation.
            for j in 1...mn {
                let kk = kn / 2
                if (kn != 2 * kk) {
                    mpmul (s2, s0, &s1, mpnw1)
                    mpeq (s1, &s2, mpnw1)
                }
                kn = kk
                if j < mn {
                    mpmul (s0, s0, &s1, mpnw1)
                    mpeq (s1, &s0, mpnw1)
                }
            }
        }
        
        ///   Compute reciprocal if N is negative.
        // 110 continue
        if n < 0 {
            mpdmc (1.0, 0, &s1, mpnw1)
            mpdiv (s1, s2, &s0, mpnw1)
            mpeq (s0, &s2, mpnw1)
        }
        
        ///   Restore original precision level.
        mproun (&s2, mpnw)
        mpeq (s2, &b, mpnw)
    } // mpnpwr

    ///   This computes the N-th root of the MPR number A and returns result in B.
    ///   N must be at least one and must not exceed 2 ^ 30.
    ///   This static func employs the following Newton-Raphson iteration, which
    ///   converges to A ^ (-1/N):
    ///   ```
    ///   X_{k+1} = X_k + (X_k / N) * (1 - A * X_k^N)
    ///   ```
    ///   The reciprocal of the final approximation to A ^ (-1/N) is the N-th root.
    ///   These iterations are performed with a maximum precision level MPNW that
    ///   is dynamically changed, approximately doubling with each iteration.
    ///   When N is large and A is very near one, the following binomial series is
    ///   employed instead of the Newton scheme:
    ///   ```
    ///   (1 + x)^(1/N)  =  1  +  x / N  +  x^2 * (1 - N) / (2/// N^2)  +  ...
    ///   ```
    ///   See the comment about the parameter NIT in MPDIVX.
    ///
    static func mpnrtr (_ a:mp_real, _ n:Int, _ b: inout mp_real, _ mpnw:Int) {
        let alt = 0.693147180559945309, cl2 = 1.4426950408889633
        let nit = 3, n30 = 1 << 30 // 2 ** 30
        var s0 = mp_init(size: mpnw+8), s1 = mp_init(size: mpnw+8), s2 = mp_init(size: mpnw+8)
        var s3 = mp_init(size: mpnw+8), f1 = mp_init(size: 9)
        
        /// End of declaration
        
        if (mpnw < 4 || a[0] < abs(a[2]) + 4 || b[0] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        let ia = sign(1, a[2])
        let na = min(Int(abs(a[2])), mpnw)
        
        if na == 0 {
//            b[1] = mpnw
//            b[2] = 0
//            b[3] = 0
//            b[4] = 0
//            b[5] = 0
            mp_zero(&b, mpnw)
            return
        }
        if ia < 0 {
            print("*** \(#function): Argument is negative.")
            mpabrt (59)
        }
        
        if n <= 0 || n > n30 {
            print("*** \(#function): Improper value of N \(n)")
            mpabrt (60)
        }
        
        ///   If N = 1 or 2, MPEQ or MPSQRT instead.
        if n == 1 {
            mpeq (a, &b, mpnw)
            return
        } else if n == 2 {
            mpsqrt (a, &b, mpnw)
            return
        }
        
        s0[0] = mpnw + 7
        s1[0] = mpnw + 7
        s2[0] = mpnw + 7
        s3[0] = mpnw + 7
        var mpnw1 = mpnw + 1
        
        ///   Set f1 = 1.
        f1[0] = 9
        f1[1] = mpnw1
        f1[2] = 1
        f1[3] = 0
        f1[4] = 1
        f1[5] = 0
        f1[6] = 0
        
        ///   Determine the least integer MQ such that 2 ^ MQ .GE. MPNW.
        var t1 = Double(mpnw)
        let mq = Int(cl2 * Double.log(t1) + 1.0 - Double(mprdfz))
        
        ///   Check how close A is to 1.
        mpsub(a, f1, &s0, mpnw1)
        if s0[2] == 0 {
            mpeq(f1, &b, mpnw)
            return
        }
        var n1 = 0
        mpmdc(s0, &t1, &n1, mpnw1)
        var n2 = Int(cl2 * Double.log(abs(t1)))
        t1 = t1 * Double.pow(0.5, Double(n2))
        n1 = n1 + n2
        
        if n1 <= -30 {
            var t2 = Double(n)
            let n2 = Int(cl2 * Double.log(t2) + 1.0 + Double(mprdfz))
            let n3 = -mpnbt * mpnw1 / n1
            if n3 < Int(1.25 * Double(n2)) {
                ///   A is so close to 1 that it is cheaper to use the binomial series.
                mpdivd (s0, t2, &s1, mpnw1)
                mpadd (f1, s1, &s2, mpnw1)
                var k = 0
                
                // 100 continue
                while true {
                    k = k + 1
                    t1 = Double(1 - k * n)
                    t2 = Double((k + 1) * n)
                    mpmuld (s1, t1, &s3, mpnw1)
                    mpdivd (s3, t2, &s1, mpnw1)
                    mpmul (s0, s1, &s3, mpnw1)
                    mpeq (s3, &s1, mpnw1)
                    mpadd (s1, s2, &s3, mpnw1)
                    mpeq (s3, &s2, mpnw1)
                    if (s1[2] != 0 && s1[3] >= -mpnw1) {
                        continue // goto 100
                    } else {
                        mpeq(s2, &s1, mpnw1)
                        mproun(&s1, mpnw)
                        mpeq(s1, &b, mpnw)
                        return
                    }
                }
            }
        }
        
        ///   Compute the initial approximation of A ^ (-1/N).
        let tn = Double(n)
        mpmdc(a, &t1, &n1, mpnw1)
        n2 = Int(Double(-n1) / tn)
        let t2 = Double.exp(-1.0 / tn * (Double.log(t1) + (Double(n1) + tn * Double(n2)) * alt))
        mpdmc(t2, n2, &s2, mpnw1)
        
        mpnw1 = 5
        var iq = 0
        
        ///   Perform the Newton-Raphson iteration described above with a dynamically
        ///   changing precision level MPNW (one greater than powers of two).
        for k in 1...mq {
            if k > 2 { mpnw1 = min(2 * mpnw1 - 2, mpnw) + 1 }
            ///  if (k > 2) mpnw1 = min (2 * mpnw1 - 1, mpnw)
            
            //110  continue
            while true {
                mpnpwr (s2, n, &s0, mpnw1)
                mpmul (a, s0, &s1, mpnw1)
                mpsub (f1, s1, &s0, mpnw1)
                mpmul (s2, s0, &s1, mpnw1)
                mpdivd (s1, Double(tn), &s0, mpnw1)
                mpadd (s2, s0, &s1, mpnw1)
                mpeq (s1, &s2, mpnw1)
                if (k == mq - nit && iq == 0) {
                    iq = 1
                } else {
                    break
                }
            }
        }
        
        ///   Take the reciprocal to give final result.
        mpdiv(f1, s2, &s1, mpnw1)
        
        ///   Restore original precision level.
        //    130 continue
        mproun(&s1, mpnw)
        mpeq(s1, &b, mpnw)
    } // mpnrtr

    ///   This returns a pseudorandom number B based the input A.
    static func mprandr (_ a:mp_real, _ b: inout mp_real, _ mpnw:Int) {
        var ia1 = mp_init(size: mpnw+7), ia2 = mp_init(size: mpnw+7), ia3 = mp_init(size: mpnw+7)
        
        if (mpnw < 4 || a[0] < abs (a[2]) + 4 || a[0] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        ia1[0] = mpnw + 7
        ia2[0] = mpnw + 7
        ia3[0] = mpnw + 7
        
        mpmuld(a, mprandx, &ia1, mpnw + 1)
        mpinfr(ia1, &ia2, &ia3, mpnw + 1)
        mpeq(ia3, &b, mpnw)
    } // mprandr
    
    ///   This performs rounding and truncation of the MPR number A.  It is called
    ///   by MPNORM, and also by other static funcs when the precision level is
    ///   modified.  It is not intended to be directly called by the user.
    ///   The parameter MPEXPMX is the absolute value of the largest exponent word
    ///   allowed for MP numbers (see system parameters at start of this module).
    static func mproun (_ a: inout mp_real, _ mpnw:Int) {
        
        if (mpnw < 4 || a[0] < abs (a[2]) + 4 || a[0] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        ///   Check for initial zeroes.
        var a2 = a[3]
        a[3] = 0
        let ia = sign(1, a[2])
        var na = min(Int(abs(a[2])), mpnw)
        let n4 = na + 4
        
        if a[4] == 0 {
            ///   Find the first nonzero word and shift the entire number left.  The length
            ///   of the result is reduced by the length of the shift.
            for i in 4...n4 {
                if a[i+1] != 0 {
                    // goto 110
                    // 110 continue
                    
                    let k = i - 3
                    for j in 3...n4-k {
                        a[j+1] = a[j+k+1]
                    }
                    
                    a2 = a2 - k
                    na = na - max (k - 2, 0)
                    if k == 2 { a[na+4] = 0 }
                    break
                } else if i == n4 {
                    a[2] = 0
                    a[3] = 0
                    a[4] = 0
                    a[5] = 0
                    return
                }
            }
        }
        
        ///   Perform rounding.
        if na == mpnw {
            if Double(a[na+4]) >= 0.5 * Double(mpbdx) { a[na+3] = a[na+3] + 1 }
            ///   Release carries as far as necessary due to rounding.
            var gotoFlag = false
            for i in stride(from: na+2, through: 3, by: -1) {
                if a[i+1] < mpbdx { gotoFlag = true; break /*goto 140*/ }
                a[i+1] = a[i+1] - mpbdx
                a[i] = a[i] + 1
            }
            
            ///   Release of carries due to rounding continued all the way to the start --
            ///   i.e. number was entirely 9"s.
            if !gotoFlag {
                a[4] = a[3]
                na = 1
                a2 = a2 + 1
            }
        }
        
        // 140 continue
        
        if a[na+3] == 0 {
            ///   At least the last mantissa word is zero.  Find the last nonzero word
            ///   and adjust the length of the result accordingly.
            var gotoFlag = false
            for i in stride(from: na+2, through: 3, by: -1) {
                if (a[i+1] != 0) { na = i - 2; gotoFlag = true; break /* goto 160 */ }
            }
            if !gotoFlag {
//                a[1] = mpnw
//                a[2] = 0
//                a[3] = 0
//                a[4] = 0
//                a[5] = 0
                mp_zero(&a, mpnw)
                return
            }
            // 160  continue
//            na = i - 2
        }
        
        ///   Check for overflow and underflow.
        if Double(a2) < -mpexpmx {
            print("*** MPROUN: Exponent underflow.")
            mpabrt (68)
        } else if Double(a2) > mpexpmx {
            print("*** MPROUN: Exponent overflow.")
            mpabrt (69)
        }
        
        ///   Check for zero.
        if a[4] == 0 {
//            a[1] = mpnw
//            a[2] = 0
//            a[3] = 0
//            a[4] = 0
//            a[5] = 0
            mp_zero(&a, mpnw)
        } else {
            a[1] = mpnw
            a[2] = sign(na, ia)
            a[3] = a2
            a[na+4] = 0
            a[na+5] = 0
        }
        
        // 170  continue
    } // mproun

    ///   This function returns 1, 0 or -1, depending on whether ra > 0, ra = 0 or ra < 0.
    static func mpsgn (_ ra:mp_real) -> Int {
        let ia = ra[2]
        if (ia == 0) {
            return 0
        } else if (ia > 0) {
            return 1
        } else {
            return -1
        }
    } // mpsgn

    ///   This computes the square root of the MPR number A and returns the result in B.
    
    ///   This function employs the following Newton-Raphson iteration, which
    ///   converges to 1 / Sqrt(A):
    
    ///    X_{k+1} = X_k + 0.5 * (1 - X_k^2 * A) * X_k
    
    ///   where the multiplication () * X_k is performed with only half of the
    ///   normal level of precision.  These iterations are performed with a
    ///   working precision level MPNW that is dynamically changed, approximately
    ///   doubling with each iteration (except that at iteration NIT before the final
    ///   iteration, the iteration is repeated without doubling the precision, in order
    ///   to enhance accuracy) .  The final iteration is performed as follows
    ///   (this is due to A. Karp):
    
    ///    Sqrt(A) = (A * X_n) + 0.5 * [A - (A * X_n)^2] * X_n  (approx.)
    
    ///   where the multiplications A * X_n and [] * X_n are performed with only
    ///   half of the final level of precision.
    static func mpsqrt (_ a:mp_real, _ b: inout mp_real, _ mpnw:Int) {
        let cl2 = 1.4426950408889633, nit = 3
        var s0 = mp_init(size: mpnw+7), s1 = mp_init(size: mpnw+7), s2 = mp_init(size: mpnw+7)
        var s3 = mp_init(size: mpnw+7)
        
        if (mpnw < 4 || a[0] < abs (a[2]) + 4 || b[0] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        let ia = sign(1, a[2])
        let na = min(Int(abs (a[2])), mpnw)
        
        if na == 0 {
//            b[1] = mpnw
//            b[2] = 0
//            b[3] = 0
//            b[4] = 0
//            b[5] = 0
            mp_zero(&b, mpnw)
            return
        }
        if ia < 0 {
            print("*** \(#function): Argument is negative.")
            mpabrt (70)
            return
        }
        
        s0[0] = mpnw + 7
        s1[0] = mpnw + 7
        s2[0] = mpnw + 7
        s3[0] = mpnw + 7
        
        ///   Determine the least integer MQ such that 2 ^ MQ .GE. MPNW.
        var t1 = Double(mpnw)
        let mq = Int(cl2 * Double.log(t1) + 1.0 - mprdfz)
        
        ///   Compute the initial approximation of 1 / Sqrt(A).
        var n = 0
        mpmdc (a, &t1, &n, mpnw)
        let n2 = -n / 2
        let t2 = Double.sqrt(t1 * Double(1 << (n + 2 * n2)))
        t1 = 1.0 / t2
        mpdmc(t1, n2, &s2, mpnw)
        mpdmc(1.0, 0, &s3, mpnw)
        
        var mpnw1 = 5
        var iq = 0
        var nw1 = mpnw1
        var nw2 = mpnw1
        
        ///   Perform the Newton-Raphson iteration described above with a dynamically
        ///   changing precision level MPNW (one greater than powers of two).
        for k in 1...mq-1 {
            if k > 2 {
                nw1 = mpnw1
                mpnw1 = min(2 * mpnw1 - 2, mpnw) + 1
                nw2 = mpnw1
            }
            
            // 100  continue
            while true {
                mpmul (s2, s2, &s0, nw2)
                mpmul (a, s0, &s1, nw2)
                mpsub (s3, s1, &s0, nw2)
                mpmul (s2, s0, &s1, nw1)
                mpmuld (s1, 0.5, &s0, nw1)
                mpadd (s2, s0, &s1, nw2)
                mpeq (s1, &s2, nw2)
                
                if k == mq - nit && iq == 0 {
                    iq = 1
                    // goto 100
                } else {
                    break
                }
            }
        }
        
        ///   Perform last iteration using Karp's trick.
        nw1 = mpnw1
        mpnw1 = min (2 * mpnw1 - 2, mpnw) + 1
        nw2 = mpnw1
        
        mpmul(a, s2, &s0, nw1)
        mpmul(s0, s0, &s1, nw2)
        mpsub(a, s1, &s3, nw2)
        mpmul(s3, s2, &s1, nw1)
        mpmuld(s1, 0.5, &s3, nw1)
        mpadd(s0, s3, &s2, nw2)
        
        ///   Restore original precision level.
        mproun(&s2, mpnw)
        mpeq(s2, &b, mpnw)
    } // mpsqrt

    ///   This routine subtracts MPR numbers A and B to yield C.
    static func mpsub (_ a:mp_real, _ b: mp_real, _ c: inout mp_real, _ mpnw:Int) {
        var s = mp_init(size: mpnw+6)
        
        if (mpnw < 4 || a[0] < abs(a[2]) + 4 || b[0] < abs(b[2]) + 4 || c[0] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        let nb = min(abs(Int(b[2])), mpnw)
        s[0] = mpnw + 6
        s[1] = mpnw
        if (b[2] == 0) {
            s[2] = 0
        } else if (b[2] > 0) {
            s[2] = -nb
        } else {
            s[2] = nb
        }
        
        s[3...nb+5] = b[3...nb+5]
//        for i = 3, nb + 5 {
//            s(i) = b(i)
//        }
        mpadd(a, s, &c, mpnw)
    } // mpsub

    // MARK: - FFT operations
    
    ///   This performs an N-point complex-to-real FFT, where N = 2^M.  X is the
    ///   double complex input array, and Y is the double precision output array.
    ///   The array X is used as a scratch array in MPFFT1, and so is overwritten.
    ///   X and Y must be dimensioned as shown below.  ISS is the sign of the FFT.
    ///   The arrays MPUU1 and MPUU2 must have been initialized by calling MPINIFFT.
    ///   This routine is not intended to be called directly by the user.
    static func mpfftcr (_ iss:Int, _ m:Int, _ n:Int, _ nsq:Int, _ x: inout [Complex], _ y: inout [Double]) {
        var dc1 = [Complex](repeating: Complex.zero, count: n/2)
        if mpuu1.isEmpty { mpinifft(mpmlxm) }
        let mx = Int(mpuu1[1].magnitude)
        
        ///   Check if input parameters are invalid.
        if ((iss != 1 && iss != -1) || m < 3 || m > mx) {
            print("*** MPFFTCR: Either the UU arrays have not been initialized",
                  "or else one of the input parameters is invalid \(iss), \(m), \(mx)")
            mpabrt (677)
        }
        
        let n1 = 1 << (m/2) // 2 ** (m / 2)
        let n2 = n / 2
        let n4 = n / 4
        let ai = Complex.i  // (0.0, 1.0)
        let half = Complex(0.5)
        
        ///   Construct the input to MPFFT1.
        dc1[1] = half * Complex((x[1] + x[n2+1]).real, (x[1] - x[n2+1]).real)
        if iss == 1 {
            dc1[n4+1] = x[n4+1].conjugate
        } else {
            dc1[n4+1] = x[n4+1]
        }
        
        let ku = n2
        if iss == 1 {
            for k in 2...n4 {
                let x1 = x[k]
                let x2 = x[n2+2-k].conjugate
                let a1 = x1 + x2
                let a2 = ai * mpuu1[k+ku] * (x1 - x2)
                dc1[k] = half * (a1 + a2)
                dc1[n2+2-k] = half * (a1 - a2).conjugate
            }
        } else {
            for k in 2...n4 {
                let x1 = x[k]
                let x2 = x[n2+2-k].conjugate
                let a1 = x1 + x2
                let a2 = ai * mpuu1[k+ku].conjugate * (x1 - x2)
                dc1[k] = half * (a1 + a2)
                dc1[n2+2-k] = half * (a1 - a2).conjugate
            }
        }
        
        ///   Perform a normal N/2-point FFT on DC1.
        // FIXME: Need to figure out this array passing
        var dc = [[Complex]](repeating: dc1, count: 1)
        var dc2 = dc
        mpfft1(iss, m - 1, n1, n2 / n1, &dc, &dc2)  // should be ... &dc1, &x)
        
        ///   Copy DC1 to Y such that DC1(k) = Y(2k-1) + i Y(2k).
        for k in 1...n/2 {
            y[2*k-1] = dc1[k].real
            y[2*k] = dc1[k].imaginary
        }
    } // mpfftcr
    
    ///   This performs an N-point real-to-complex FFT, where N = 2^M.  X is the
    ///   double precision input array, and Y is the double complex output array.
    ///   The arrays MPUU1 and MPUU2 must have been initialized by calling MPINIFFT.
    ///   This routine is not intended to be called directly by the user.
    static func mpfftrc (_ iss:Int, _ m:Int, _ n:Int, _ nsq:Int, _ x:[Double], _ y: inout [Complex]) {
        
        //    implicit none
        //    integer, intent(in):: is, m, n, nsq
        //    integer k, ku, mx, n1, n2, n4
        //    real (mprknd), intent(in):: x(n)
        //    complex (mprknd), intent(out):: y(n/2+nsq*mpnsp1+1)
        //    complex (mprknd) dc1(n/2), ai, a1, a2, z1, z2
        var dc1 = [Complex](repeating: Complex.zero, count: n/2)
        if mpuu1.isEmpty { mpinifft(mpmlxm) }
        let mx = Int(mpuu1[1].magnitude)
        
        ///   Check if input parameters are invalid.
        if ((iss != 1 && iss != -1) || m < 3 || m > mx) {
            print("*** MPFFTRC: either the UU arrays have not been initialized",
                  "or else one of the input parameters is invalid \(iss), \(m), \(mx)")
            mpabrt (677)
        }
        
        let n1 = 1 << (m / 2) // 2 ** (m / 2)
        let n2 = n / 2
        let n4 = n / 4
        let ai = Complex(0.0, -1.0)
        
        ///   Copy X to DC1 such that DC1(k) = X(2k-1) + i X(2k).
        for k in 1...n2 {
            dc1[k] = Complex(x[2*k-1], x[2*k])
        }
        
        ///   Perform a normal N/2-point FFT on DC1.
        // FIXME: Need to figure out this array passing
        var dc = [[Complex]](repeating: dc1, count: 1)
        var dc2 = dc
        mpfft1(iss, m - 1, n1, n2 / n1, &dc, &dc2)  // was ...  &dc1, &y)
        
        ///   Reconstruct the FFT of X.
        y[1] = Complex(2.0 * (dc1[1].real + dc1[1].imaginary), 0.0)
        let two = Complex(2.0)
        if iss == 1 {
            y[n4+1] = two * dc1[n4+1]
        } else {
            y[n4+1] = two * dc1[n4+1].conjugate
        }
        y[n2+1] = Complex(2.0 * (dc1[1].real - dc1[1].imaginary), 0.0)
        let ku = n2
        
        if iss == 1 {
            for k in 2...n4 {
                let z1 = dc1[k]
                let z2 = dc1[n2+2-k].conjugate
                let a1 = z1 + z2
                let a2 = ai * mpuu1[k+ku] * (z1 - z2)
                y[k] = a1 + a2
                y[n2+2-k] = (a1 - a2).conjugate
            }
        } else {
            for k in 2...n4 {
                let z1 = dc1[k]
                let z2 = dc1[n2+2-k].conjugate
                let a1 = z1 + z2
                let a2 = ai * mpuu1[k+ku].conjugate * (z1 - z2)
                y[k] = a1 + a2
                y[n2+2-k] = (a1 - a2).conjugate
            }
        }
        
    } // mpfftrc

    ///   This routine performs a complex-to-complex FFT.  IS is the sign of the
    ///   transform, N = 2^M is the size of the transform.  N1 = 2^M1 and N2 = 2^M2,
    ///   where M1 and M2 are defined as below.  X is the input and output array,
    ///   and Y is a scratch array.  X must have at N, and Y at least N + N1\*MPNSP1,
    ///   double complex cells.  The arrays MPUU1 and MPUU2 must have been
    ///   initialized by calling MPINIFFT.  This routine is not intended to be called
    ///   directly by the user.
    ///   This employs the two-pass variant of the "four-step" FFT.  See the
    ///   article by David H. Bailey in J. of Supercomputing, March 1990, p. 23-35.
    static func mpfft1 (_ iss:Int, _ m:Int, _ n1:Int, _ n2:Int, _ x: inout [[Complex]], _ y: inout [[Complex]]) {
        
        //    implicit none
        //    integer, intent(in):: is, m, n1, n2
        //    integer i, iu, j, j2, k, ku, m1, m2, nr1, nr2
        //    complex (mprknd), intent(inout):: x(n1,n2)
        //    complex (mprknd), intent(out):: y(n2+mpnsp1,n1)
        //    complex (mprknd) z1(mpnrow+mpnsp1,n1), z2(mpnrow+mpnsp1,n1)
        
        let m1 = (m + 1) / 2
        let m2 = m - m1
        let nr1 = min (n1, mpnrow)
        let nr2 = min (n2, mpnrow)
        let ku = Int(mpuu2[m].magnitude)
        
        let zi = [Complex](repeating: Complex.zero, count: mpnrow+mpnsp1)
        var z1 = [[Complex]](repeating: zi, count: n1)
        var z2 = z1
        
        for i in stride(from: 0, to: n1, by: nr1) {
            ///   Copy NR1 rows of X (treated as a N1 x N2 complex array) into Z1.
            for j in 1...n2 {
                for k in 1...nr1 {
                    z1[k][j] = x[i+k][j]
                }
            }
            
            ///   Perform NR1 FFTs, each of length N2.
            mpfft2(iss, nr1, m2, n2, &z1, &z2)
            
            ///   Multiply the resulting NR1 x N2 complex block by roots of unity and
            ///   store transposed into the appropriate section of Y.
            let iu = i + ku - n1 - 1
            if (iss == 1) {
                for j in 1...n2 {
                    for k in 1...nr1 {
                        y[j][i+k] = mpuu2[iu+k+j*n1] * z1[k][j]
                    }
                }
            } else {
                for j in 1...n2 {
                    for k in 1...nr1 {
                        y[j][i+k] = mpuu2[iu+k+j*n1].conjugate * z1[k][j]
                    }
                }
            }
        }
        
        for i in stride(from: 0, to: n2, by: nr2) {
            ///   Copy NR2 rows of the Y array into Z2.
            for j in 1...n1 {
                for k in 1...nr2 {
                    z2[k][j] = y[i+k][j]
                }
            }
            
            ///   Perform NR2 FFTs, each of length N1.
            mpfft2 (iss, nr2, m1, n1, &z2, &z1)
            
            ///   Copy NR2 x N1 complex block back into X array.  It's a little more
            ///   complicated if M is odd.
            if m.isMultiple(of: 2) {
                for j in 1...n1 {
                    for k in 1...nr2 {
                        x[i+k][j] = z2[k][j]
                    }
                }
            } else {
                for j in 1...n1/2 {
                    let j2 = 2 * j - 1
                    for k in 1...nr2 {
                        x[i+k][j] = z2[k][j2]
                        x[i+k+n2][j] = z2[k][j2+1]
                    }
                }
            }
        }
    } // mpfft1

    ///   This performs NS simultaneous N-point complex-to-complex FFTs, where
    ///   N = 2^M.  X is the input and output array, and Y is a scratch array.
    ///   The arrays MPUU1 and MPUU2 must have been initialized by calling MPINIFFT.
    ///   This routine is not intended to be called directly by the user.
    static func mpfft2 (_ iss:Int, _ ns:Int, _ m:Int, _ n:Int, _ x: inout [[Complex]], _ y: inout [[Complex]]) {
        ///   Perform the second variant of the Stockham FFT.
        for l in stride(from: 1, through: m, by: 2) {
            mpfft3(iss, l, ns, m, n, &x, &y)
            if l == m { // goto 100
                // 100 continue
                ///   Copy Y to X.
                for j in 1...n {
                    for i in 1...ns {
                        x[i][j] = y[i][j]
                    }
                }
                return
            }
            mpfft3(iss, l + 1, ns, m, n, &y, &x)
        }
    } // mpfft2

    ///   This performs the L-th iteration of the second variant of the Stockham FFT
    ///   on the NS vectors in X.  X is input/output, and Y is a scratch array.
    ///   The arrays MPUU1 and MPUU2 must have been initialized by calling MPINIFFT.
    ///   This routine is not intended to be called directly by the user.
    static func mpfft3 (_ iss:Int, _ l:Int, _ ns:Int, _ m:Int, _ n:Int, _ x: inout [[Complex]], _ y: inout [[Complex]]) {
        ///   Set initial parameters.
        if mpuu1.isEmpty { mpinifft(mpmlxm) }
        let n1 = n / 2
        let lk = 1 << (l-1) // 2 ** (l - 1)
        let li = 1 << (m-1) // 2 ** (m - l)
        let lj = 2 * lk
        let ku = li + 1
        var u1: Complex
        for i in 0..<li {
            let i11 = i * lk + 1
            let i12 = i11 + n1
            let i21 = i * lj + 1
            let i22 = i21 + lk
            if iss == 1 {
                u1 = mpuu1[i+ku]
            } else {
                u1 = mpuu1[i+ku].conjugate
            }
            
            for k in 0..<lk {
                for j in 1...ns {
                    let x1 = x[j][i11+k]
                    let x2 = x[j][i12+k]
                    y[j][i21+k] = x1 + x2
                    y[j][i22+k] = u1 * (x1 - x2)
                }
            }
        }
    } // mpfft3

    ///   This computes the root of unity arrays UU1 and UU2, which are required by
    ///   the FFT routines, and places this data in the proper arrays defined in
    ///   module MPFUNA.  MPNW is the largest precision level (in words) that will be
    ///   subsequently used for this run.
    static func mpinifft (_ mpnw:Int) {
        let cl2 = 1.4426950408889633
        
        ///  Determine sizes for FFT arrays.  Three words are added to mpnw, since many
        ///  routines in MPFUND in particular increase the working precision upon entry.
        let nwds = mpnw + 3
        let d1 = 2.0 * Double(nwds + 1)
        let m = Int(cl2 * Double.log(d1) + 1.0 - mprdfz)
        let mq = m + 2
        let nq = 1 << mq // 2 ** mq
        
        if mq + nq > mplfftx {
            print("*** MPINIFFT: Insufficient space for arrays mpuu1 and mpuu2.",
                  "At least \(mq+nq) double complex cells must be allocated for each of",
                  "these arrays in module mpfuna. See documentation for details.")
            mpabrt (91)
        }
        
        mpuu1[1] = Complex(mq)
        var ku = 2
        var ln = 1
        let pi = Double.pi
        
        for _ in 1...mq {
            let t1 = pi / Double(ln)
            for i in 0..<ln {
                let ti = Double(i) * t1
                mpuu1[i+ku] = Complex(Double.cos(ti), Double.sin(ti))
            }
            
            ku = ku + ln
            ln = 2 * ln
        }
        
        /// write (6, 2) ku - 1
        /// 2 format ("MPINIFFT: Size of table mpuu1 =",i10)
        
        ku = mq + 1
        mpuu2[1] = Complex(mq)
        
        for k in 2...mq {
            mpuu2[k] = Complex.zero
        }
        
        for k in 2..<mq {
            mpuu2[k] = Complex(ku)
            let mm = k
            let nn = 1 << mm // 2 ** mm
            let mm1 = (mm + 1) / 2
            let mm2 = mm - mm1
            let nn1 = 1 << mm1 // 2 ** mm1
            let nn2 = 1 << mm2 // 2 ** mm2
            let tpn = 2.0 * pi / Double(nn)
            
            for j in 0..<nn2 {
                for i in 0..<nn1 {
                    let iu = ku + i + j * nn1
                    let t1 = tpn * Double(i * j)
                    mpuu2[iu] = Complex(Double.cos(t1), Double.sin(t1))
                }
            }
            
            ku = ku + nn
        }
        /// write (6, 3) ku - 1
        /// 3 format ("MPINIFFT: Size of table mpuu2 =",i10)
    } // mpinifft
    
    ///   This computes the linear convolution of A and B, returning the result
    ///   in C.  If IQ is 1, then it is presumed B = A; if IQ = 2, then A != B.
    ///   NSQ is a spacing parameter, which should be set to more than sqrt (3*n).
    static func mplconv (_ iq:Int, _ n:Int, _ nsq:Int, _ a:[Double], _ b:[Double], _ c: inout [Double]) {
        var d1 = [Double](repeating: 0, count: 8*n+2)
        var d2 = [Double](repeating: 0, count: 8*n+2)
        var d3 = [Double](repeating: 0, count: 8*n+2)
        var dc1 = [Complex](repeating: Complex.zero, count: 4*n+nsq*mpnsp1+3)
        var dc2 = [Complex](repeating: Complex.zero, count: 4*n+nsq*mpnsp1+3)
        let cl2 = 1.4426950408889633, ffterrmx = 0.375
        
        var t1 = Double(n)
        let m1 = Int(cl2 * Double.log(t1) + 1.0 - mprdfz)
        let n1 = 1 << m1 // 2 ** m1
        let m2 = m1 + 1
        let n2 = 2 * n1
        let n4 = 2 * n2
        let nm = min (2 * n, n2)
        
        if abs(iq) == 1 {
            ///   Compute the square of A -- only one forward FFT is needed.
            for i in 1...n {
                d1[i] = a[i]
            }
            
            for i in n+1...n2 {
                d1[i] = 0
            }
            
            ///   Perform a forward real-to-complex FFT on the vector in A.
            mpfftrc (1, m2, n2, nsq, d1, &dc1)
            
            ///   Square the resulting complex vector.
            for i in 1...n1+1 {
                dc1[i] *= dc1[i]
            }
        } else {
            ///   Compute the product of A and B -- two forward FFTs are needed.
            for i in 1...n {
                d1[i] = a[i]
                d2[i] = b[i]
            }
            
            for i in n+1...n2 {
                d1[i] = 0
                d2[i] = 0
            }
            
            ///   Perform forward real-to-complex FFTs on the vectors in A and B.
            mpfftrc (1, m2, n2, nsq, d1, &dc1)
            mpfftrc (1, m2, n2, nsq, d2, &dc2)
            
            ///   Multiply the resulting complex vectors.
            for i in 1...n1+1 {
                dc1[i] = dc1[i] * dc2[i]
            }
        }
        
        ///   Perform an inverse complex-to-real FFT on the resulting data.
        mpfftcr(-1, m2, n2, nsq, &dc1, &d3)
        
        ///   Divide by N4 and round to nearest whole number.
        let an = 1.0 / Double(n4)
        var c0 = 0.0
        
        for i in 1...nm {
            t1 = an * Double(d3[i])
            let t2 = t1.rounded()
            c[i] = t2
            c0 = max(c0, abs(t2 - t1))
        }
        
        if c0 > ffterrmx {
            print("*** MPLCONV: excessive rounding error = \(c0)")
            mpabrt (55)
        }
    } // mplconv
    
    ///   This routine multiplies MP numbers A and B to yield the MP product C,
    ///   using a FFT-convolution technique.  Before calling MPMULX, the arrays
    ///   UU1 and UU2 must be initialized by calling MPINIFFT.  For modest levels
    ///   of precision, use MPMUL.
    static func mpmulx (_ a:mp_real, _ b:mp_real, _ c: inout mp_real, _ mpnw:Int) {
        var d1 = [Double](repeating: 0, count: 4*mpnw+21)
        var d2 = [Double](repeating: 0, count: 4*mpnw+21)
        var d3 = [Double](repeating: 0, count: 4*mpnw+41)
        var d = mp_init(size: mpnw+9)
        
        if (mpnw < 4 || a[0] < abs (a[2]) + 4 || b[0] < abs (b[2]) + 4 || c[0] < mpnw + 6) {
            print("*** MPMULX: uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        let ia = sign(1, a[2])
        let ib = sign(1, b[2])
        let na = min (Int(abs (a[2])), mpnw)
        let nb = min (Int(abs (b[2])), mpnw)
        let nc = min (na + nb, mpnw)
        let nn = 4 * max (na, nb)
        let nx = Int(Double.sqrt(4.0 * Double(nn)) + mprdfz)
        
        ///   Divide each word of A into four 15-bit chunks.
        for i in 0..<na {
            var i1 = a[i+4]
            var i2 = i1 >> 45
            d1[4*i] = Double(i2)
            i1 = i1 - (i2 << 45)
            i2 = i1 >> 30
            d1[4*i+1] = Double(i2)
            i1 = i1 - (i2 << 30)
            i2 = i1 >> 15
            d1[4*i+2] = Double(i2)
            i1 = i1 - (i2 << 15)
            d1[4*i+3] = Double(i1)
        }
        
        for i in 4*na..<nn {
            d1[i] = 0
        }
        
        ///   Divide each word of B into four 15-bit chunks.
        for i in 0..<nb {
            var i1 = b[i+4]
            var i2 = i1 >> 45
            d2[4*i] = Double(i2)
            i1 = i1 - (i2 << 45)
            i2 = i1 >> 30
            d2[4*i+1] = Double(i2)
            i1 = i1 - (i2 << 30)
            i2 = i1 >> 15
            d2[4*i+2] = Double(i2)
            i1 = i1 - (i2 << 15)
            d2[4*i+3] = Double(i1)
        }
        
        for i in 4*nb..<nn {
            d2[i] = 0
        }
        
        ///   Perform linear convolution.
        mplconv(2, nn, nx, d1, d2, &d3)
        
        ///   Release carries.
        var i0 = 0
        for i in stride(from: min(4 * nc + 16, 2 * nn - 1), through: 0, by: -1) {
            i0 = i0 + Int(d3[i])
            let i1 = i0 >> 15
            let i2 = i0 - (i1 << 15)
            d1[i] = Double(i2)
            i0 = i1
        }
        
        ///  Recombine words, with proper offset.
        d[0] = 0
        d[1] = 0
        d[2] = 0
        d[3] = 0
        d[4] = (i0 << 45) + (Int(d1[0]) << 30) + (Int(d1[1]) << 15) + Int(d1[2])
        
        for i in 1...nc+3 {
            d[i+4] = (Int(d1[4*i-1]) << 45) + (Int(d1[4*i]) << 30) + (Int(d1[4*i+1]) << 15) + Int(d1[4*i+2])
        }
        
        d[0] = mpnw + 6
        d[1] = mpnw
        d[2] = sign(nc, ia * ib)
        d[3] = a[3] + b[3] + 1
        
        ///   Fix up the result.
        // d1[0] = mpnw + 6   // why is this being done?  It is never used
        mpnorm(&d, &c, mpnw)
    } // mpmulx

    // MARK: - String to MPReal conversions
    
    ///  Converts the numeric string in *a* into the MPR number *b*.
    ///  Restrictions: (a) no embedded blanks; (b) a leading digit (possibly
    ///  zero) must be present; and (c) a period must be present.  An exponent
    ///  (with "d" or "e") may optionally follow the numeric value.
    static func mpctomp (_ a:String, _ b: inout mp_real, _ mpnw: Int) {
        let a = a.trimmingCharacters(in: .whitespacesAndNewlines)
        let lexpmx = 9
//        let digits = "0123456789"
        let d10w = Double.pow(10.0, Double(mpndpw))
        var f = mp_init(size: 9)
        var s0 = mp_init(size: mpnw+7), s1 = mp_init(size: mpnw+7), s2 = mp_init(size: mpnw+7)
        
        /// write (6, *) "mpctomp: a, n, mpnw =", n, mpnw
        /// write (6, "(100a1)") "X",(a(i),i=1,n),"X"
        
        func abort(_ err:Int) {
            print("*** \(#function): Syntax error in input string; code = \(err)",
                  "Restrictions: (a) no embedded blanks; (b) a leading digit (possibly",
                  "zero) must be present; and (c) a period must be present.  An exponent",
                  "(with \"d\" or \"e\") may optionally follow the numeric value.")
            mpabrt(41)
        }
        
        if (mpnw < 4 || b[0] < mpnw + 6) {
            // write (mpldb, 1)
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        s0[0] = mpnw + 7
        s1[0] = mpnw + 7
        s2[0] = mpnw + 7
        f[0] = 9
        f[1] = mpnw
        
        for i in 2...8 { f[i] = 0 }
        
        let mpnw1 = mpnw + 1
        var kde = -1
        let kend = a.count-1
        var kexpend = -1
        var kexpst = -1
        var kexpsgn = -1
        var knumend1 = -1
        var knumend2 = -1
        var knumst1 = -1
        var knumst2 = -1
        var kper = -1
        var ksgn = -1
        let kstart = 0
        
        ///   Locate:
        ///     kstart = index of first nonblank character.
        ///     kend = index of last nonblank character.
        if a.isEmpty {
            ///   Input is completely blank.
            abort(1)
        }
        
        ///   Scan input for:
        ///     kde = index of "d" or "e".
        ///     kexpend = index of end of exponent.
        ///     kexpst = index of start of exponent.
        ///     kespsgn = index of sign of exponent.
        ///     knumend1 = index of end of numeric part prior to period.
        ///     knumend2 = index of end of numeric part after period.
        ///     knumst1 = index of start of numeric part prior to period.
        ///     knumst2 = index of start of numeric part after period.
        ///     kper = index of period.
        ///     ksgn = index of sign of number.
        
        for (i,ch) in a.enumerated() {
            if ch == " " {
                abort(2)
            } else if ["+", "-"].contains(ch) {
                if i == kstart {
                    ksgn = i
                } else if kde > 0 && kexpsgn < 0 && kexpst < 0 && i < kend {
                    kexpsgn = i
                } else {
                    abort(3)
                }
            } else if ["E", "D"].contains(ch.uppercased()) {
                if kde < 0 && kper >= 0 && i < kend {
                    kde = i
                    knumend2 = i - 1
                } else {
                    abort(4)
                }
            } else if ch == "." {
                if kper < 0 && kde < 0 && knumst1 >= 0 && knumst2 < 0 {
                    kper = i
                    knumend1 = i - 1
                } else {
                    abort(5)
                }
            } else if ch.isWholeNumber {
                if knumst1 < 0 {
                    knumst1 = i
                } else if kper >= 0 && knumst2 < 0 && kde < 0 {
                    knumst2 = i
                } else if kde >= 0 && kexpst < 0 {
                    kexpst = i
                }
                if i == kend {
                    if knumst2 >= 0 && kde < 0 {
                        knumend2 = i
                    } else if kexpst >= 0 {
                        kexpend = i
                    } else {
                        abort(6)
                    }
                }
            } else {
                abort(7)
            }
        }
        
        /// write (6, *) "kde, kend, kexpend, kexpst =", kde, kend, kexpend, kexpst
        /// write (6, *) "kexpsgn, numend1, knumend2, knumst1 =", kexpsgn, knumend1, knumend2, knumst1
        /// write (6, *) "knumst2, kper, ksgn, kstart =", knumst2, kper, ksgn, kstart
        
        ///   Decode exponent.
        var iexp = 0.0
        let chs = Array(a)
        if kexpst > 0 {
            let lexp = kexpend - kexpst + 1
            if lexp > lexpmx {
                print("*** \(#function): exponent string is too long.")
                mpabrt (85)
            }
            let ca = String(chs[kexpst...kexpst+lexp-1])
            iexp = mpdigin (ca, lexp)
            if kexpsgn > 0 {
                if chs[kexpsgn] == "-" { iexp = -iexp }
            }
        }
        
        ///   Determine lengths of two sections of number.
        let lnum1 = knumend1 - knumst1 + 1
        let lnum2:Int
        if knumst2 > 0 {
            lnum2 = knumend2 - knumst2 + 1
        } else {
            lnum2 = 0
        }
        let lnum = lnum1 + lnum2
        
        /// write (6, *) "iexp, lnum1, lnum2 =", iexp, lnum1, lnum2
        
        ///   Determine the number of chunks of digits and the left-over.
        let n1 = lnum / mpndpw
        let n2 = lnum % mpndpw
        
        ///   Construct first (left-over) portion, right-justified in CA.
        var ca = ""
        var ix = knumst1-1
        for _ in 1...n2 {
            ix = ix + 1
            if ix == kper { ix = ix + 1 }
            ca.append(chs[ix])
        }
        
        var t1 = mpdigin(ca, mpndpw)
        if t1 > 0 {
            f[2] = 1
            f[3] = 0
            f[4] = Int(t1)
        } else {
            f[2] = 0
            f[3] = 0
            f[4] = 0
        }
        mpeq(f, &s0, mpnw1)
        
        ///   Process remaining chunks of digits.
        for _ in 1...n1 {
            ca = ""
            for _ in 1...mpndpw {
                ix = ix + 1
                if ix == kper { ix = ix + 1 }
                ca.append(chs[ix])
            }
            
            t1 = mpdigin(ca, mpndpw)
            if t1 > 0 {
                f[2] = 1
                f[3] = 0
                f[4] = Int(t1)
            } else {
                f[2] = 0
                f[3] = 0
                f[4] = 0
            }
            mpmuld (s0, d10w, &s1, mpnw1)
            mpadd (s1, f, &s0, mpnw1)
        }
        
        ///  Correct exponent.
        iexp = iexp - Double(lnum2)
        f[2] = 1
        f[3] = 0
        f[4] = 10
        mpnpwr(f, Int(iexp), &s1, mpnw1)
        mpmul(s0, s1, &s2, mpnw1)
        if (ksgn > 0) {
            if (chs[ksgn] == "-") { s2[2] = -s2[2] }
        }
        
        ///   Restore original precision and exit.
        mproun(&s2, mpnw)
        mpeq(s2, &b, mpnw)
        
        /// write (6, *) "mpctomp: output ="
        /// mpout (6, 420, 400, b, mpnw)
        
    } // static func mpctomp

    ///   This converts the string CA of nonblank length N to double precision.
    ///   CA may only be modest length and may only contain digits.  Blanks are ignored.
    ///   This is intended for internal use only.
    static func mpdigin (_ ca:String, _ n:Int) -> Double {
        let digits = Array("0123456789")
        var d1 = 0.0
        for ch in ca {
            if ch != " " {
                if let k = digits.firstIndex(of: ch) {
                    d1 = 10.0 * d1 + Double(k)
                } else {
                    print("*** \(#function): non-digit in character string = \(ch)")
                    mpabrt (86)
                }
            }
        }
        
        return d1
    } // mpdigin
    
    static func aint(_ x: Double) -> Double { x.rounded(.towardZero) }
    
    ///   This converts the double precision input A to a character(32) string of
    ///   nonblank length N.  A must be a whole number, and N must be sufficient
    ///   to hold it.  This is intended for internal use only.
    static func mpdigout (_ a:Double, _ n:Int) -> String {
//        let digits = Array("0123456789")
//        var ca = ""
        let d1 = Int(abs(a))
        let ca = String(d1)
//        for _ in stride(from: n, through: 1, by: -1) {
//            let d2 = aint(d1 / 10.0)
//            let k = Int(d1 - 10.0 * d2)
//            d1 = d2
//            ca.insert(digits[k], at: ca.startIndex)
//        }
        return "".padding(toLength: n-ca.count, withPad: "0", startingAt: 0) + ca
    } // mpdigout
    
    static func sign(_ v:Int, _ i:Int) -> Int { Int(v.magnitude) * i.signum() }
    static func sign(_ v:Double, _ i:Double) -> Double { i.sign == .minus ? -v.magnitude : v.magnitude }

    // MARK: - MPReal to string conversions
    
    ///   Converts the MPR number A into character form in the character[1] array B.
    ///   NB (input) is the length of the output string, and ND (input) is the
    ///   number of digits after the decimal point.  The format is analogous to
    ///   Fortran E format.  The result is left-justified among the NB cells of B.
    ///   The condition NB >= ND + 10 must hold or an error message will result.
    ///   NB cells must be available in array B.
    static func mpeformat (_ a:mp_real, _ nb:Int, _ nd:Int, _ b: inout String, _ mpnw:Int) {
        let digits = Array("0123456789")
        let d10w = Double.pow(10.0, Double(mpndpw))
        var f = mp_init(size: 9)
        var s0 = mp_init(size: mpnw+7), s1 = mp_init(size: mpnw+7), s2 = mp_init(size: mpnw+7)
        
        if (mpnw < 4 || a[0] < abs(a[2]) + 4 || nb < nd + 10) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        let ia = sign(1, a[2])
        let na = min(abs(a[2]), mpnw)
        s0[0] = mpnw + 7
        s1[0] = mpnw + 7
        s2[0] = mpnw + 7
        let mpnw1 = mpnw + 1
        
        ///   Set f = 10.
        f[0] = 9
        f[1] = mpnw1
        f[2] = 1
        f[3] = 0
        f[4] = 10
        f[5] = 0
        f[6] = 0
        
        ///   Determine power of ten for exponent, and scale input to within 1 and 10.
        var nexp:Int
        if na > 0 {
            var aa = Double(a[4])
            if na >= 2 { aa = aa + Double(a[5]) / Double(mpbdx) }
            let t1 = Double.log10(2.0) * Double(mpnbt * a[3]) + Double.log10(aa)
            if t1 >= 0.0 {
                nexp = Int(t1)
            } else {
                nexp = Int(t1 - 1.0)
            }
            
            if nexp == 0 {
                mpeq(a, &s1, mpnw1)
            } else if nexp > 0 {
                mpnpwr(f, nexp, &s0, mpnw1)
                mpdiv(a, s0, &s1, mpnw1)
            } else if nexp < 0 {
                mpnpwr (f, -nexp, &s0, mpnw1)
                mpmul (a, s0, &s1, mpnw1)
            }
            
            ///   If we didn't quite get it exactly right, multiply or divide by 10 to fix.
            while true {
                if s1[3] < 0 {
                    nexp = nexp - 1
                    mpmuld (s1, 10.0, &s0, mpnw1)
                    mpeq (s0, &s1, mpnw1)
                } else if (s1[4] >= 10) {
                    nexp = nexp + 1
                    mpdivd(s1, 10.0, &s0, mpnw1)
                    mpeq(s0, &s1, mpnw1)
                } else {
                    break
                }
            }
            s1[2] = abs(s1[2])
        } else {
            nexp = 0
            mpeq(a, &s1, mpnw1)
        }
        
        ///   Insert sign and first digit.
        var b2 = ""
        if ia == -1 {
            b2 += "-"
        }
        var an:Double
        if na > 0 {
            an = Double(s1[4])
        } else {
            an = 0.0
        }
        var ca = mpdigout(an, 1)
        b2 += ca
        b2 += "."
        let ixp = b2.count
        
        ///   Set f = an.
        f[0] = 9
        f[1] = mpnw1
        f[2] = 1
        f[3] = 0
        f[4] = Int(an)
        f[5] = 0
        f[6] = 0
        mpsub(s1, f, &s0, mpnw1)
        mpmuld(s0, d10w, &s1, mpnw1)
        
        ///   Calculate the number of remaining chunks.
        let nl = nd / mpndpw + 1
        
        ///   Insert the digits of the remaining words.
        for _ in 1...nl {
            if (s1[2] != 0 && s1[3] == 0) {
                an = Double(s1[4])
                f[2] = 1
                f[3] = 0
                f[4] = Int(an)
            } else {
                f[2] = 0
                f[3] = 0
                f[4] = 0
                an = 0.0
            }
            
            let ca = mpdigout(an, mpndpw)
            b2 += ca
            
            mpsub(s1, f, &s0, mpnw1)
            mpmuld(s0, d10w, &s1, mpnw1)
        }
        
        ///   Round the result.
        if b2.count >= nd+1 {
            var breakFlag = false
            var b2c = Array(b2)
            let d = b2c[nd+1]
        //    let i1 = index (digits, b2(nd+1)) - 1
            if var i1 = digits.firstIndex(of: d), i1 >= 5 {
                ///   Perform rounding, beginning at the last digit (position ND).  If the rounded
                ///   digit is 9, set to 0, } repeat at position one digit to left.  Continue
                ///   rounding if necessary until the decimal point is reached.
                for i in stride(from: nd, through: ixp+1, by: -1) {
                    let d = b2c[i]
                    let i2 = digits.firstIndex(of: d)!
                    if i2 <= 8 {
                        b2.append(digits[i2+1])
                        breakFlag = true
                        break
                      //   goto 180
                    } else {
                        b2 += "0"
                    }
                }
                
                ///   We have rounded up all digits to the right of the decimal point.  If the
                ///   digit to the left of the decimal point is a 9, then set that digit to 1
                ///   and increase the exponent by one; otherwise increase that digit by one.
                if !breakFlag {
                    if b2c[ixp-1] == "9" {
                        b2c[ixp-1] = "1"
                        nexp = nexp + 1
                    } else {
                        i1 = digits.firstIndex(of: b2c[ixp-1])!
                        b2c[ixp-1] = digits[i1+1]
                    }
                    b2 = String(b2c)  // convert back to string
                }
            }
        }
        
        // 180 continue
        
        ///   Done with mantissa.  Insert exponent.
        b2 += "e"
        if nexp < 0 {
            b2 += "-"
        }
        ca = mpdigout(Double(abs(nexp)), 10)
        
        // drop leading zeros
        while ca.hasPrefix("0") && ca.count > 1 {
            ca.removeFirst()
        }
        
        b2 += ca
        while b2.count < nb {
            b2 += " "
        }
        
        ///   Copy entire b2 array to B.
        b = b2
    } // mpeformat

    ///   Converts the MPR number A into character form in the character[1] array B.
    ///   NB (input) is the length of the output string, and ND (input) is the
    ///   number of digits after the decimal point.  The format is analogous to
    ///   Fortran F format; the result is right-justified among the NB cells of B.
    ///   The condition NB >= ND + 10 must hold or an error message will result.
    static func mpfformat (_ a:mp_real, _ nb:Int, _ nd:Int, _ b: inout String, _ mpnw:Int) {
        if (mpnw < 4 || a[0] < abs(a[2]) + 4 || nb < nd + 10) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        var nb2 = nb + 20
        var b2 = ""
        mpeformat(a, nb2, nb, &b2, mpnw+1)
        
        ///   Trim off trailing blanks.
        while b2.hasSuffix(" ") {
            b2.removeLast()
        }
        nb2 = b2.count
        
        ///   Look for the "e" in B2.
        let b2c = Array(b2) // convert to Character array
        var k = b2c.firstIndex(of: "e") ?? -1
        if k == -1 {
            print("*** \(#function): Syntax error in output of mpeformat")
            mpabrt (84)
        }
        
        ///   Check the sign of the exponent.
        let ixp:Int
        k = k + 1
        if b2c[k] == "-" {
            ixp = -1
            k = k + 1
        } else {
            ixp = 1
        }
        
        ///   Copy the exponent into CA.
        var j = 0
        var ca = ""
        for i in k..<b2c.count {
            j = j + 1
            if j <= 16 { ca.append(b2c[i]) }
        }
        
        let t1 = mpdigin(ca, j)
        
        ///   Check if there is enough space in the output array for all digits.
        ///   Not necessary with Swift
//        if (t1 + nd + 3 > nb) {
//            for i = 1, nb {
//                b(i) = "*"
//            }
//            return
//        }
        let nexp = ixp * Int(t1)
        
        ///   Insert the sign of the number, if any.
        b = ""
        var ch = b2.removeFirst()
        if ch == "-" {
            b += "-"
            ch = b2.removeFirst()
        }
        
        if nexp == 0 {
            ///   Exponent is zero.  Copy first digit, period and ND more digits.
            for _ in 1...nd+2 {
                b.append(ch)
                ch = b2.removeFirst()
            }
        } else if nexp > 0 {
            ///   Exponent is positive.  Copy first digit, skip the period, then copy
            ///   nexp digits.
            b.append(ch)
            ch = b2.removeFirst() // drop the period
            for _ in 1...nexp {
                ch = b2.removeFirst()
                b.append(ch)
            }
            
            ///   Insert the period.
            b += "."
            
            ///   Copy nd more digits.
            for _ in 1...nd {
                if b2.isEmpty { break }
                ch = b2.removeFirst()
                b.append(ch)
            }
        } else {
            ///   Exponent is negative.  Insert a zero, then a period, then nexp - 1
            ///   zeroes, then the first digit, then the remaining digits up to ND total
            ///   fractional digits.
            b += "0."
            for _ in stride(from: 1, to: nexp, by: 1) {
                b += "0"
            }
            
            b.append(ch)
            ch = b2.removeFirst() // drop the period
            
            for _ in nexp...nd {
                ch = b2.removeFirst()
                b.append(ch)
            }
        }
        
        ///   Right-justify in field.
        while b.count < nd {
            b = " " + b
        }
    } // mpfformat
    
    // MARK: - Text file reading operation
    
    ///   This routine reads the MPR number A from a string *ifile*.  The *ifile*
    ///   string is consumed during the formation of *a*; however, only those
    ///   characters that are actually required for *a*.  The digits of *a*
    ///   may span more than one line, provided that a "\" appears at the end of
    ///   a line to be continued (any characters after the "\" on the same line
    ///   are ignored).  Embedded blanks are allowed anywhere.
    ///   An exponent with "e" or "d" may optionally follow the numeric value.
    ///   Note: Normally the *ifile* is read from a URL in one operation before
    ///   this function is invoked.
    static func mpinp (_ ifile: inout String, _ a: inout mp_real, _ mpnw:Int) {
        
        func getLine() -> String {
            ifile = ifile.trimmingCharacters(in: .whitespacesAndNewlines)
            let lines = ifile.components(separatedBy: .newlines)
            while lines.count > 0, let line = lines.first {
                ifile = String(ifile.dropFirst(line.count))  // remove line from input
                let trimmedLine = line.trimmingCharacters(in: .whitespaces)
                if !trimmedLine.isEmpty { return trimmedLine }
            }
            print("*** \(#function): End-of-file encountered.")
            mpabrt (72)
            return ""
        }
        
        let validc = Array(" 0123456789+-.dDeE")
        if (mpnw < 4 || a[0] < mpnw + 6) {
            print("*** \(#function): uninitialized or inadequately sized arrays")
            mpabrt (99)
        }
        
        let lncx = mpnw * Int(mpdpw + 1) + 1000
        var line1 = getLine()
        
        ///   Scan input line(s), looking for valid characters.
        var chr1 = ""
        while !line1.isEmpty {
            let ch = line1.removeFirst()
            if ch == "\\" { line1 = getLine(); continue }
            let i1 = validc.firstIndex(of: ch) ?? -1
            if i1 == -1 && ch != " " {
                print("*** MPINP: Invalid input character = \(ch)")
                mpabrt (87)
            } else if ch != " " {
                if chr1.count < lncx {
                    chr1.append(ch)
                }
            }
        }
        
        mpctomp(chr1, &a, mpnw)
    } // static func mpinp

}


// MARK: - Some constants used to define the mp numbers

extension MPReal {
    
    /// Set the default standard and medium precision levels (in digits) here.
    static public var mpipl = 2500
    static public var mpiplm = 250
    
    /// Do not change the following code (in normal usage).
    static public let mpwds = Int(Double(mpipl) / mpdpw + 2.0)
    static public let mpwds6 = mpwds + 6
    static public let mpwdsm = Int(Double(mpiplm) / mpdpw + 2.0)
    static public let mpwdsm6 = mpwdsm + 6
    
    /// Precision level in words for usage of mpmulx.
    /// Default 1111 words corresponds to approximately 20,000 digits.
    /// If higher precision is required, mpinifft or mpinit.
    static public let mpmlxm = 1111
    
    /// Precision in words of precomputed log(2), pi and egamma below.
    static public let mpl2pi = 1111
    
    /// Largest n such that 10^n <= 2^53. See usage in MPFUNC.
    static public let mpndpw = 15
    
    /// Number of mantissa bits per long integer word in MP numbers.
    static public let mpnbt = 60
    
    /// Maximum length of certain input character strings. See usage in mpinp of module MPFUNC.
    static public let mpnstr = 2048
    
    /// Length of output lines. See usage in mpout of MPFUNC.
    static public let mpoutl = 80
    
    /// No quad-precision Float available
    static public let mprknd2 = -1
    
    /// 2^mpnbt, the radix for MP numbers.
    static public let mpbdx = 1 << mpnbt
    
    /// DP approximation to number of digits per mantissa word.
    static public let mpdpw = 18.061799739838871713
    
    /// Largest permissible exponent, corresponding to a maximum binary exponent of 2^31.
    static public let mpexpmx = Double.pow(2.0, 31.0) / Double(mpnbt)
    
    /// "Fuzz" for comparing DP values.
    static let mprdfz = Double.pow(0.5, 50.0)
    
    /// Multiplier for pseudorandom number generator.
    static let mprandx = 5501758857861179.0
    
    /// Max size of log(2) and pi arrays. This is large enough to
    /// support precision levels over one million digits.
    static public let mplconx = 100000
    
    /// Max size of FFT arrays. This is large enough to support
    /// precision levels over one million digits.
    static public let mplfftx = 2100000
    
    /// Row & spacing parameters for FFT.
    static let mpnrow = 16
    static let mpnsp1 = 2
    static let mpnsp2 = 9
}
