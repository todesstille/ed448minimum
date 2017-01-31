package ed448

import . "gopkg.in/check.v1"

var (
	primeSerial = [fieldBytes]byte{
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xfe, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
	}

	one  = [fieldBytes]byte{1}
	zero = [fieldBytes]byte{}
)

func (s *Ed448Suite) Test_PointMul(c *C) {
	resultZero := PointMul(zero, testValue)
	resultValueTimes1 := PointMul(one, testValue)

	c.Assert(resultZero, DeepEquals, zero[:])
	c.Assert(resultValueTimes1, DeepEquals, testValue[:])

	val := PointMul(primeSerial, one)
	c.Assert(val, IsNil)

	val = PointMul(one, primeSerial)
	c.Assert(val, IsNil)

	val = PointMul(primeSerial, primeSerial)
	c.Assert(val, IsNil)
}

func (s *Ed448Suite) Test_PointAdd(c *C) {
	valuePlusOne := [fieldBytes]byte{
		0x04, 0x44, 0x58, 0xab, 0x92, 0xc2, 0x78,
		0x23, 0x55, 0x8f, 0xc5, 0x8d, 0x32, 0xc2,
		0x6c, 0x21, 0x90, 0x36, 0xd6, 0xae, 0x49,
		0xdb, 0x4e, 0xc4, 0xe9, 0x23, 0xca, 0x7c,
		0xff, 0xff, 0xff, 0x1f, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0x2f, 0xff, 0xff, 0xff, 0xff, 0xff,
		0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x3f,
	}

	resultAddZero := PointAddition(zero, testValue)
	resultAddOne := PointAddition(one, testValue)

	c.Assert(resultAddZero, DeepEquals, testValue[:])
	c.Assert(resultAddOne, DeepEquals, valuePlusOne[:])

	val := PointAddition(primeSerial, primeSerial)
	c.Assert(val, IsNil)

	val = PointAddition(primeSerial, zero)
	c.Assert(val, IsNil)

	val = PointAddition(zero, primeSerial)
	c.Assert(val, IsNil)
}

func (s *Ed448Suite) Test_ScalarCopyEquals(c *C) {
	a := NewScalar([fieldBytes]byte{})
	b := a.Copy()
	c.Assert(a.Equals(b), Equals, true)
}
