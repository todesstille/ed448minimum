package ed448

type twPNiels struct {
	n *twNiels
	z *bigNumber
}

func (p *twPNiels) toExtendedPoint() *twExtendedPoint {
	eu := &bigNumber{}
	q := &twExtendedPoint{
		&bigNumber{},
		&bigNumber{},
		&bigNumber{},
		&bigNumber{},
	}
	eu.add(p.n.b, p.n.a)
	q.y.sub(p.n.b, p.n.a)
	q.t.mul(q.y, eu)
	q.x.mul(p.z, q.y)
	q.y.mul(p.z, eu)
	q.z.square(p.z)

	return q
}

type twNiels struct {
	a, b, c *bigNumber
}

func newNielsPoint(a, b, c [fieldBytes]byte) *twNiels {
	return &twNiels{
		a: mustDeserialize(serialized(a)),
		b: mustDeserialize(serialized(b)),
		c: mustDeserialize(serialized(c)),
	}
}

func (p *twNiels) conditionalNegate(neg word) {
	p.a.conditionalSwap(p.b, neg)
	p.c = p.c.conditionalNegate(neg)
}
