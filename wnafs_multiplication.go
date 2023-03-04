package ed448

import (
	"math/bits"
)

type smvtControl struct {
	power, addend int
}

const bOver16 = uint(2)

func recodeWNAF2(control []smvtControl, s *scalar, tableBits uint) word {
	tableSize := 446/(tableBits+1) + 3
	position := tableSize - 1

	control[position].power = -1
	control[position].addend = 0
	position--

	current := uint64(s[0] & 0xFFFF)
	mask := uint32((1 << (tableBits + 1)) - 1)

	for w := uint(1); w < 30; w++ {
		if w < 28 {
			// Refill the 16 high bits of current
			current += uint64(uint32((s[w/bOver16] >> uint((16 * (w % bOver16)))) << 16))
		}

		for current&0xFFFF != 0 {
			pos := uint32(bits.TrailingZeros32(uint32(current)))
			odd := uint32(current) >> pos
			delta := int32(odd & mask)
			if odd&(1<<(tableBits+1)) != 0 {
				delta -= (1 << (tableBits + 1))
			}
			current = uint64(int64(current) - int64(delta<<pos))
			control[position].power = int(pos) + 16*(int(w)-1)
			control[position].addend = int(delta)
			position--
		}
		current >>= 16
	}

	position++
	n := uint(tableSize - position)
	for i := uint(0); i < n; i++ {
		control[i] = control[i+position]
	}
	return word(n - 1)
}

func (p *twExtendedPoint) prepareFixedWindow(nTable int) []*twPNiels {
	pOriginal := p.copy()
	pn := p.copy().doubleInternal(false).toPNiels()
	out := make([]*twPNiels, nTable)
	out[0] = pOriginal.toPNiels()
	for i := 1; i < nTable; i++ {
		pOriginal.addProjectiveNielsToExtended(pn, false)
		out[i] = pOriginal.toPNiels()
	}
	return out[:]
}

func decafPrepareWNAFTable(dst []*twPNiels, p *twExtendedPoint, tableSize uint) {
	dst[0] = p.toPNiels()

	if tableSize == 0 {
		return
	}

	p.doubleInternal(false)

	twOp := p.toPNiels()

	p.addProjectiveNielsToExtended(dst[0], false)
	dst[1] = p.toPNiels()

	for i := 2; i < 1<<tableSize; i++ {
		p.addProjectiveNielsToExtended(twOp, false)
		dst[i] = p.toPNiels()
	}
}

func decafDoubleNonSecretScalarMul(p *twExtendedPoint, scalarPre, scalarVar *scalar) *twExtendedPoint {
	tableBitsVar := uint(3) // DECAF_WNAF_VAR_TABLE_BITS
	tableBitsPre := uint(5) // DECAF_WNAF_FIXED_TABLE_BITS

	var controlVar [115]smvtControl // nbitsVar/(tableBitsVar+1)+3
	var controlPre [77]smvtControl  // nbitsPre/(tableBitsPre+1)+3

	recodeWNAF2(controlPre[:], scalarPre, tableBitsPre)
	recodeWNAF2(controlVar[:], scalarVar, tableBitsVar)

	var precmpVar [8]*twPNiels
	decafPrepareWNAFTable(precmpVar[:], p, tableBitsVar)

	contp := 0
	contv := 0

	index := controlVar[0].addend >> 1

	i := controlVar[0].power
	out := &twExtendedPoint{
		&bigNumber{0x00},
		&bigNumber{0x00},
		&bigNumber{0x00},
		&bigNumber{0x00},
	}

	if i > controlPre[0].power {
		out = precmpVar[index].toExtendedPoint()
		contv++
	} else if i == controlPre[0].power && i >= 0 {
		out = precmpVar[index].toExtendedPoint()
		out.addNielsToExtended(decafWnafsTable[controlPre[0].addend>>1], i != 0)
		contv++
		contp++
	} else {
		i = controlPre[0].power
		out = decafWnafsTable[controlPre[0].addend>>1].toExtended()
		contp++
	}

	if i < 0 {
		out.setIdentity()
		return out
	}

	for i--; i >= 0; i-- {

		cv := i == controlVar[contv].power
		cp := i == controlPre[contp].power

		out.doubleInternal(i != 0 && !(cv || cp))

		if cv {
			if controlVar[contv].addend > 0 {
				out.addProjectiveNielsToExtended(precmpVar[controlVar[contv].addend>>1], (i != 0 && !cp))
			} else {
				out.subProjectiveNielsFromExtendedPoint(precmpVar[(-controlVar[contv].addend)>>1], (i != 0 && !cp))
			}
			contv++
		}

		if cp {
			if controlPre[contp].addend > 0 {
				out.addNielsToExtended(decafWnafsTable[controlPre[contp].addend>>1], i != 0)
			} else {
				out.subNielsFromExtendedPoint(decafWnafsTable[(-controlPre[contp].addend)>>1], i != 0)
			}
			contp++
		}
	}
	return out
}
