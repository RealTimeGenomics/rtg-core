/*
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.util.arithcode;

/**
 * <P>Interface for an adaptive statistical model of a stream to be used
 * as a basis for arithmetic coding and decoding. As in {@link
 * java.io.InputStream}, bytes are coded as integers in the range
 * <code>0</code> to <code>255</code> and <code>EOF</code> is provided
 * as a constant and coded as <code>-1</code>.  In addition,
 * arithmetic coding requires an integer <code>ESCAPE</code> to code
 * information about the model structure.
 *
 * <P> During encoding, a series of calls will be made to
 * <code>escaped(symbol)</code> where <code>symbol</code> is a byte
 * encoded as an integer in the range 0 to 255 or <code>EOF</code>,
 * and if the result is <code>true</code>, a call to
 * <code>interval(ESCAPE)</code> will be made and the process repeated
 * until a call to <code>escaped(symbol)</code> returns
 * <code>false</code>, at which point a call to
 * <code>interval(symbol)</code> is made and the underlying model is
 * updated.
 *
 * <P> During decoding, a call to <code>total()</code> will be made
 * and then a call to <code>pointToSymbol(count)</code>.  If the
 * result is <code>ESCAPE</code>, the process is repeated.  If the
 * result is a byte encoded as an integer in the range <code>0</code>
 * to <code>255</code> or <code>EOF</code>, the symbol is returned and
 * the underlying model is updated.
 *
 * <P>The probability model required for arithmetic coding is
 * cumulative.  For each outcome, rather than returning a probability,
 * an interval is provided to the coder.  As is usual for arithmetic
 * coding, an interval in <code>[0,1]</code> is represented by three
 * integers, where a low count, a high count, and a total count pick
 * out the interval <code>[low/total,high/total)</code>.
 *
 * <P> For more details, see <a href="../../../tutorial.html">The
 * Arithmetic Coding Tutorial</a>.
 *
 * @version 1.1
 */
public interface ArithCodeModel {

  /**
   * Use the model to encode one symbol into the encoder.
   * @param encoder used to do encoding.
   * @param symbol to be encoded.
   */
  void encode(ArithEncoder encoder, int symbol);

  /**
   * Use the model to decode a single symbol from the decoder.
   * @param decoder used to do decoding.
   * @return the symbol.
   */
  int decode(ArithDecoder decoder);
}
