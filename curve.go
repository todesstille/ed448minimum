package ed448

// Deserializes an array of bytes (little-endian) into an array of words
func bytesToWords(dst []word, src []byte) {
	wordBytes := uint(wordBits / 8)
	srcLen := uint(len(src))

	dstLen := uint((srcLen + wordBytes - 1) / wordBytes)
	if dstLen < uint(len(dst)) {
		panic("wrong dst size")
	}

	for i := uint(0); i*wordBytes < srcLen; i++ {
		out := word(0)
		for j := uint(0); j < wordBytes && wordBytes*i+j < srcLen; j++ {
			out |= word(src[wordBytes*i+j]) << (8 * j)
		}

		dst[i] = out
	}
}
