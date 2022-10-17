# taken from https://github.com/quantumgizmos/bias_tailored_qldpc
import ldpc.protograph as pt
import numpy as np
from bposd.css import css_code
from bposd.stab import stab_code


def Identity(n):
    return pt.identity(n)


class lifted_hgp(css_code):
    def __init__(self, lift_parameter, a, b=None):

        """
        Generates the lifted hypergraph product of the protographs a and b
        """
        self.a = a

        self.a_m, self.a_n = self.a.shape

        if b is None:
            self.b = pt.copy(self.a)
        else:
            self.b = b

        self.b_m, self.b_n = self.b.shape

        self.hx1_proto = np.kron(self.a, Identity(self.b_n))
        self.hx2_proto = np.kron(Identity(self.a_m), self.b.T)
        self.hx_proto = pt.hstack([self.hx1_proto, self.hx2_proto])

        self.hz1_proto = np.kron(Identity(self.a_n), self.b)
        self.hz2_proto = np.kron(self.a.T, Identity(self.b_m))
        self.hz_proto = pt.hstack([self.hz1_proto, self.hz2_proto])

        self.lift_parameter = lift_parameter

        super().__init__(self.hx_proto.to_binary(lift_parameter), self.hz_proto.to_binary(lift_parameter))

    @property
    def protograph(self):
        px = pt.vstack([pt.zeros(self.hz_proto.shape), self.hx_proto])
        pz = pt.vstack([self.hz_proto, pt.zeros(self.hx_proto.shape)])
        return pt.hstack([px, pz])

    @property
    def hx1(self):
        return self.hx1_proto.to_binary(self.lift_parameter)

    @property
    def hx2(self):
        return self.hx2_proto.to_binary(self.lift_parameter)

    @property
    def hz1(self):
        return self.hz1_proto.to_binary(self.lift_parameter)

    @property
    def hz2(self):
        return self.hz2_proto.to_binary(self.lift_parameter)


class bias_tailored_lifted_product(stab_code):
    def __init__(self, lift_parameter, a, b=None):
        lhgp = lifted_hgp(lift_parameter, a, b)

        # Hadamard rotation
        temp1 = pt.hstack([pt.zeros(lhgp.hx1_proto.shape), lhgp.hz2_proto])
        temp2 = pt.hstack([lhgp.hx1_proto, pt.zeros(lhgp.hz2_proto.shape)])
        self.hx_proto = pt.vstack([temp1, temp2])
        temp1 = pt.hstack([lhgp.hz1_proto, pt.zeros(lhgp.hx2_proto.shape)])
        temp2 = pt.hstack([pt.zeros(lhgp.hz1_proto.shape), lhgp.hx2_proto])
        self.hz_proto = pt.vstack([temp1, temp2])

        super().__init__(self.hx_proto.to_binary(lift_parameter), self.hz_proto.to_binary(lift_parameter))

    @property
    def protograph(self):
        px = pt.vstack([pt.zeros(self.hz_proto.shape), self.hx_proto])
        pz = pt.vstack([self.hz_proto, pt.zeros(self.hx_proto.shape)])
        return pt.hstack([px, pz])
