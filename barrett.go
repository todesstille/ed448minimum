package ed448

type barrettPrime struct {
	wordsInP word
	pShift   word
	lowWords []word
}

var curvePrimeOrder = barrettPrime{
	wordsInP: wordsInP,
	pShift:   pShift,
	lowWords: lowWords,
}

func barrettDeserializeAndReduce(dst []word, serial []byte, p *barrettPrime) {
	wordLen := wordBits / 8
	size := (len(serial) + wordLen - 1) / wordLen
	if size < int(p.wordsInP) {
		size = int(p.wordsInP)
	}

	tmp := make([]word, size)
	bytesToWords(tmp[:], serial[:])
	barrettReduce(tmp[:], 0, p)

	for i := word(0); i < p.wordsInP; i++ {
		dst[i] = tmp[i]
	}
}

func barrettReduce(dst []word, carry word, p *barrettPrime) {
	for wordsLeft := word(len(dst)); wordsLeft >= p.wordsInP; wordsLeft-- {
		//TODO PERF unroll
		for repeat := 0; repeat < 2; repeat++ {
			mand := dst[wordsLeft-1] >> p.pShift
			dst[wordsLeft-1] &= (word(1) << p.pShift) - 1

			if p.pShift != 0 && repeat == 0 {
				if wordsLeft < word(len(dst)) {
					mand |= dst[wordsLeft] << (wordBits - p.pShift)
					dst[wordsLeft] = 0
				} else {
					mand |= carry << (wordBits - p.pShift)
				}
			}

			carry = widemac(
				dst[wordsLeft-p.wordsInP:wordsLeft],
				p.lowWords, mand, 0)
		}
	}

	cout := addExtPacked(dst, dst[:p.wordsInP], p.lowWords, lmask)

	if p.pShift != 0 {
		cout = (cout << (wordBits - p.pShift)) + (dst[p.wordsInP-1] >> p.pShift)
		dst[p.wordsInP-1] &= word(1)<<p.pShift - 1
	}

	// mask = carry-1: if no carry then do sub, otherwise don't
	subExtPacked(dst, dst[:p.wordsInP], p.lowWords, cout-1)
}

func addExtPacked(dst, x, y []word, mask word) word {
	carry := sdword(0)
	for i := 0; i < len(y); i++ {
		carry += sdword(x[i]) + sdword(y[i]&mask)
		dst[i] = word(carry)
		carry >>= wordBits
	}

	for i := len(y); i < len(x); i++ {
		carry += sdword(x[i])
		dst[i] = word(carry)
		carry >>= wordBits
	}

	return word(carry)
}

func subExtPacked(dst, x, y []word, mask word) word {
	carry := sdword(0x00)
	for i := 0; i < len(y); i++ {
		carry += sdword(x[i]) - (sdword(y[i]) & sdword(mask))
		dst[i] = word(carry)
		carry >>= wordBits
	}

	for i := len(y); i < len(x); i++ {
		carry += sdword(x[i])
		dst[i] = word(carry)
		carry >>= wordBits
	}

	return word(carry)
}

// TODO Is this the same as mulAddVWW_g() ?
func widemac(accum []word, mier []word, mand, carry word) word {
	for i := 0; i < len(mier); i++ {
		product := dword(mand) * dword(mier[i])
		product += dword(accum[i])
		product += dword(carry)

		accum[i] = word(product)
		carry = word(product >> wordBits)
	}

	for i := len(mier); i < len(accum); i++ {
		sum := dword(carry) + dword(accum[i])
		accum[i] = word(sum)
		carry = word(sum >> wordBits)
	}

	return carry
}
