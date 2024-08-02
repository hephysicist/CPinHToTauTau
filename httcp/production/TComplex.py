import os

from columnflow.util import maybe_import
np = maybe_import("numpy")
ak = maybe_import("awkward")


class TComplex:
    def __init__(self, re=0, im=0, polar=False):
        """
           Mimicing ROOT TComplex class
           this re or im can be any numbers / arrays ( awkward, numpy )
        """

        if isinstance(re, ak.Array):
            if isinstance(im, (int, float)):
                im = im * ak.ones_like(re)
        elif isinstance(re, np.ndarray):
            if isinstance(im, (int, float)):
                im = im * np.ones_like(re)
                
        self.set_complex_numbers(re, im, polar)

        
    def set_complex_numbers(self, re, im, polar):        
        if polar:
            if isinstance(re, (int, float)):
                if re < 0:
                    print("Warning: Modulo of a complex number should be >= 0, taking the absolute value.")
                    re = np.abs(re)
                self.fRe = re * np.cos(im)
                self.fIm = re * np.sin(im)
            elif isinstance(re, ak.Array):
                self.fRe = ak.where(re < 0, np.abs(re) * np.cos(im), re * np.cos(im))
                self.fIm = ak.where(re < 0, np.abs(re) * np.sin(im), re * np.sin(im))
            elif isinstance(re, np.ndarray):
                self.fRe = np.where(re < 0, np.abs(re) * np.cos(im), re * np.cos(im))
                self.fIm = np.where(re < 0, np.abs(re) * np.sin(im), re * np.sin(im))
        else:
            self.fRe = re
            self.fIm = im

    def Re(self):
        return self.fRe

    def Im(self):
        return self.fIm

    def Rho(self):
        return np.sqrt(self.fRe**2 + self.fIm**2)

    def Rho2(self):
        return self.fRe**2 + self.fIm**2

    def Theta(self):
        if isinstance(self.fRe, (int, float)):
            return np.arctan2(self.fIm, self.fRe) if (self.fIm or self.fRe) else 0

        elif isinstance(self.fRe, np.ndarray):
            dummy = np.zeros_like(self.fRe)
            return np.where(mask, np.arctan2(self.fIm, self.fRe), dummy)
        elif isinstance(self.fRe, ak.Array):
            dummy = ak.zeros_like(self.fRe)
            return ak.where(mask, np.arctan2(self.fIm, self.fRe), dummy)
        else:
            raise RuntimeError("no specific type")

    def Sqrt(self):
        rho = np.sqrt(self.Rho())
        theta = 0.5 * self.Theta()
        return TComplex(rho * np.cos(theta), rho * np.sin(theta), polar=True)

    def Exp(self):
        return TComplex(np.exp(self.fRe), self.fIm, polar=True)

    def Log(self):
        return TComplex(0.5 * np.log(self.Rho2()), self.Theta())

    def Sin(self):
        return TComplex(np.sin(self.fRe) * np.cosh(self.fIm), np.cos(self.fRe) * np.sinh(self.fIm), polar=True)

    def Cos(self):
        return TComplex(np.cos(self.fRe) * np.cosh(self.fIm), -np.sin(self.fRe) * np.sinh(self.fIm), polar=True)

    def Tan(self):
        cc = self.Cos()
        return self.Sin() * cc.Conjugate() / cc.Rho2()

    def ASin(self):
        return -TComplex.I() * TComplex.Log(TComplex.I() * self + np.sign(1., self.Im()) * TComplex.Sqrt(1. - self**2))

    def ACos(self):
        return -TComplex.I() * TComplex.Log(self + np.sign(1., self.Im()) * TComplex.Sqrt(self**2 - 1.))

    def ATan(self):
        return -0.5 * TComplex.I() * TComplex.Log((1. + TComplex.I() * self) / (1. - TComplex.I() * self))

    def SinH(self):
        return TComplex(np.sinh(self.fRe) * np.cos(self.fIm), np.cosh(self.fRe) * np.sin(self.fIm), polar=True)

    def CosH(self):
        return TComplex(np.cosh(self.fRe) * np.cos(self.fIm), np.sinh(self.fRe) * np.sin(self.fIm), polar=True)

    def TanH(self):
        cc = self.CosH()
        return self.SinH() * cc.Conjugate() / cc.Rho2()

    def ASinH(self):
        return TComplex.Log(self + np.sign(1., self.Im()) * TComplex.Sqrt(self**2 + 1.))

    def ACosH(self):
        return TComplex.Log(self + np.sign(1., self.Im()) * TComplex.Sqrt(self**2 - 1.))

    def ATanH(self):
        return 0.5 * TComplex.Log((1. + self) / (1. - self))

    def Conjugate(self):
        return TComplex(self.fRe, -self.fIm)

    def __mul__(self, other):
        if isinstance(other, TComplex):
            return TComplex(self.fRe * other.fRe - self.fIm * other.fIm, self.fRe * other.fIm + self.fIm * other.fRe)
        elif isinstance(other, (int, float, np.ndarray, ak.Array)):
            return TComplex(self.fRe * other, self.fIm * other)
        else:
            raise TypeError(f"Unsupported type for multiplication: {type(other)}")

    def __rmul__(self, other):
        return self.__mul__(other)

    
    def __add__(self, other):
        if isinstance(other, TComplex):
            return TComplex(self.fRe + other.fRe, self.fIm + other.fIm)
        elif isinstance(other, (int, float, np.ndarray, ak.Array)):
            return TComplex(self.fRe + other, self.fIm)
        else:
            raise TypeError(f"Unsupported type for addition: {type(other)}")

    def __radd__(self, other):
        return self.__add__(other)

        
    def __truediv__(self, other):
        if isinstance(other, TComplex):
            rho2 = other.Rho2()
            return TComplex((self.fRe * other.fRe + self.fIm * other.fIm) / rho2,
                            (-self.fRe * other.fIm + self.fIm * other.fRe) / rho2)
        elif isinstance(other, (int, float, np.ndarray, ak.Array)):
            return TComplex(self.fRe / other, self.fIm / other)
        else:
            raise TypeError(f"Unsupported type for division: {type(other)}")

    def __sub__(self, other):
        if isinstance(other, TComplex):
            return TComplex(self.fRe - other.fRe, self.fIm - other.fIm)
        elif isinstance(other, (int, float, np.ndarray, ak.Array)):
            return TComplex(self.fRe - other, self.fIm)
        else:
            raise TypeError(f"Unsupported type for subtraction: {type(other)}")

    def __neg__(self):
        return TComplex(-self.fRe, -self.fIm)

    def __pos__(self):
        return self

    def __pow__(self, other):
        if isinstance(other, TComplex):
            lrho = np.log(self.Rho())
            theta = self.Theta()
            return TComplex(np.exp(lrho * other.fRe - theta * other.fIm), lrho * other.fIm + theta * other)

    def __iadd__(self, other):
        if isinstance(other, TComplex):
            self.fRe += other.fRe
            self.fIm += other.fIm
        elif isinstance(other, (int, float, np.array, ak.Array)):
            self.fRe += other
        else:
            raise TypeError(f"Unsupported type for addition: {type(other)}")
        return self
