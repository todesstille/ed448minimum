package ed448

import (
	"crypto/rand"
	"crypto/sha512"
)

// Curve is the interface that wraps the basic curve methods.
// TODO It would be better with the use of privateKey and publicKey types.
type Curve interface {
	GenerateKeys() (priv [privKeyBytes]byte, pub [pubKeyBytes]byte, ok bool)
	Sign(priv [privKeyBytes]byte, message []byte) (signature [signatureBytes]byte, ok bool)
	Verify(signature [signatureBytes]byte, message []byte, pub [pubKeyBytes]byte) (valid bool)
	ComputeSecret(private [privKeyBytes]byte, public [pubKeyBytes]byte) (secret [sha512.Size]byte)
}

type curveT struct{}

var (
	curve = &curveT{}
)

type DecafCurve interface {
	GenerateKeys() (priv [privKeyBytes]byte, pub [pubKeyBytes]byte, ok bool)
	Sign(priv [privKeyBytes]byte, message []byte) (signature [signatureBytes]byte, ok bool)
	Verify(signature [signatureBytes]byte, message []byte, pub [pubKeyBytes]byte) (valid bool, err error)
}

type decafCurveT struct{}

var (
	decafCurve = &decafCurveT{}
)

// NewDecafCurve returns a Curve.
func NewDecafCurve() DecafCurve {
	return decafCurve
}

// GenerateKeys generates a private key and its correspondent public key.
func (ed *decafCurveT) GenerateKeys() (priv [privKeyBytes]byte, pub [pubKeyBytes]byte, ok bool) {
	privKey, err := ed.decafGenerateKeys(rand.Reader)
	ok = err == nil

	copy(priv[:], privKey[:])
	copy(pub[:], privKey.publicKey())

	return
}

// Signs a message using the provided private key and returns the signature.
func (ed *decafCurveT) Sign(priv [privKeyBytes]byte, message []byte) (signature [signatureBytes]byte, ok bool) {
	pk := privateKey(priv)

	signature, err := ed.decafSign(message, &pk)
	ok = err == nil
	return
}

// Verify a signature does correspond a message by a public key.
func (ed *decafCurveT) Verify(signature [signatureBytes]byte, message []byte, pub [pubKeyBytes]byte) (valid bool, err error) {
	pk := publicKey(pub)

	valid, err = ed.decafVerify(signature, message, &pk)
	if err != nil {
		return false, err
	}
	return valid, nil
}
