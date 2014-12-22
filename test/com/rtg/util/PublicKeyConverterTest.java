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
package com.rtg.util;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.security.KeyPairGenerator;
import java.security.interfaces.RSAPublicKey;
import java.util.Arrays;

import junit.framework.TestCase;

/**
 * Tests corresponding class
 */
public class PublicKeyConverterTest extends TestCase {

  public void testToFrom() throws Exception {
    final KeyPairGenerator kpg = KeyPairGenerator.getInstance("RSA");
    kpg.initialize(512);
    final RSAPublicKey rsapub = (RSAPublicKey) kpg.generateKeyPair().getPublic();
    final ByteArrayOutputStream bkey = new ByteArrayOutputStream();
    PublicKeyConverter.fromPublicKey(rsapub, bkey);
    final ByteArrayInputStream bkeyIn = new ByteArrayInputStream(bkey.toByteArray());
    RSAPublicKey result = PublicKeyConverter.toPublicKey(bkeyIn);
    assertEquals(rsapub, result);
    assertEquals(rsapub.getPublicExponent(), result.getPublicExponent());
    assertEquals(rsapub.getModulus(), result.getModulus());
  }

  private static final byte[] BYTE_ARRAYS_SRC = {0, 0, 0, 8, //size1
                                                0, 0, 3, 0, 5, 0, 7, 8, //value1
                                                0, 0, 0, 4, //size2
                                                1, 0, 3, 4}; //value2
  private static final byte[][] BYTE_ARRAYS_RESULT = {{3, 0, 5, 0, 7, 8},
                                                      {1, 0, 3, 4}};

  public void testToByteArrays() throws IOException {
    byte[][] result = PublicKeyConverter.toByteArrays(new ByteArrayInputStream(BYTE_ARRAYS_SRC));
    assertEquals(BYTE_ARRAYS_RESULT.length, result.length);
    assertTrue(Arrays.equals(BYTE_ARRAYS_RESULT[0], result[0]));
    assertTrue(Arrays.equals(BYTE_ARRAYS_RESULT[1], result[1]));
  }
}
