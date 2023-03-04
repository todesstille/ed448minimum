package ed448

func mask(a, b *bigNumber, mask word) {
	a[0] = word(mask) & b[0]
	a[1] = word(mask) & b[1]
	a[2] = word(mask) & b[2]
	a[3] = word(mask) & b[3]
	a[4] = word(mask) & b[4]
	a[5] = word(mask) & b[5]
	a[6] = word(mask) & b[6]
	a[7] = word(mask) & b[7]
	a[8] = word(mask) & b[8]
	a[9] = word(mask) & b[9]
	a[10] = word(mask) & b[10]
	a[11] = word(mask) & b[11]
	a[12] = word(mask) & b[12]
	a[13] = word(mask) & b[13]
	a[14] = word(mask) & b[14]
	a[15] = word(mask) & b[15]
}
