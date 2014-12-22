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

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.OutputStream;
import java.math.BigInteger;
import java.security.KeyFactory;
import java.security.NoSuchAlgorithmException;
import java.security.interfaces.RSAPublicKey;
import java.security.spec.InvalidKeySpecException;
import java.security.spec.RSAPublicKeySpec;

import com.rtg.util.gzip.GzipUtils;

/**
 * Turns our internal public key format into a Java <code>PublicKey</code>
 * and vice versa
 */
public final class PublicKeyConverter {

  private PublicKeyConverter() {
  }

  /**
   * Converts an <code>RSAPublicKey</code> to our internal representation which is binary:
   * <code>[int32-modulus-size][big-endian-modulus][int32-exponent-size][big-endian-exponent]</code>
   * @param key key to convert
   * @param dest destination for data
   * @throws IOException If an IO error occurred
   */
  public static void fromPublicKey(final RSAPublicKey key, final OutputStream dest) throws IOException {
    final DataOutputStream dos = new DataOutputStream(dest);
    try {
      final byte[] modBytes = key.getModulus().toByteArray();
      dos.writeInt(modBytes.length);
      dos.write(modBytes);
      final byte[] expBytes = key.getPublicExponent().toByteArray();
      dos.writeInt(expBytes.length);
      dos.write(expBytes);
    } finally {
      dos.flush();
    }
  }

  //.NET tizz over these
  private static byte[] stripLeadingZeros(final byte[] b) {
    int start = 0;
    for (int i = 0; i < b.length; i++) {
      if (b[i] == 0) {
        start = i + 1;
      } else {
        break;
      }
    }
    final byte[] ret = new byte[b.length - start];
    System.arraycopy(b, start, ret, 0, ret.length);
    return ret;
  }

  /**
   * Takes our internal RSA public key format and turns it into a <code>RSAPublicKey</code>
   * @param key stream containing our internal format
   * @return public key object
   * @throws IOException if an IO Error occurs
   */
  public static RSAPublicKey toPublicKey(final InputStream key) throws IOException {
    final DataInputStream dis = new DataInputStream(key);
    final int modLength = dis.readInt();
    final byte[] modBytes = new byte[modLength];
    dis.readFully(modBytes);
    final int expLength = dis.readInt();
    final byte[] expBytes = new byte[expLength];
    dis.readFully(expBytes);

    final BigInteger modulus = new BigInteger(modBytes);
    final BigInteger publicExponent = new BigInteger(expBytes);

    final KeyFactory factory;
    try {
      factory = KeyFactory.getInstance("RSA");
    } catch  (final NoSuchAlgorithmException nsae) {
      throw new RuntimeException("Couldn't find RSA algorithm", nsae);
    }
    final RSAPublicKeySpec keyspec = new RSAPublicKeySpec(modulus, publicExponent);
    try {
      return (RSAPublicKey) factory.generatePublic(keyspec);
    } catch (final InvalidKeySpecException ikse) {
      throw new IllegalArgumentException("Provided RSA values invalid.", ikse);
    }
  }

  /**
   * Takes our internal RSA public key format and turns it into byte arrays
   * to apply to a .NET <code>RSAParameters</code> object
   * Note: Not used by our Java code.
   * @param key input data
   * @return modulus and exponent byte arrays respectively
   * @throws IOException If an IO Error occurs
   */
  public static byte[][] toByteArrays(final InputStream key) throws IOException {
    final DataInputStream dis = new DataInputStream(key);
    final int modLength = dis.readInt();
    final byte[] modBytes = new byte[modLength];
    dis.readFully(modBytes);
    final int expLength = dis.readInt();
    final byte[] expBytes = new byte[expLength];
    dis.readFully(expBytes);
    return new byte[][] {stripLeadingZeros(modBytes), stripLeadingZeros(expBytes)};
  }

  /**
   * Internal use only, to ease the production of our internal public key
   * from a serialized one.
   * @param args first value is filename of gzipped serialized RSA public key
   * @throws java.io.IOException If it happens
   * @throws java.lang.ClassNotFoundException very unlikely
   */
  public static void main(final String[] args) throws IOException, ClassNotFoundException {
    if (args.length < 2) {
      System.err.println("Usage: [public-key-file] [output-file]");
      return;
    }
    final File pubFile = new File(args[0]);
    if (!pubFile.exists()) {
      System.err.println("public key file does not exist");
      return;
    }
    final File outFile = new File(args[1]);
    if (outFile.exists()) {
      System.err.println("Output file already exists");
      return;
    }
    final RSAPublicKey key;
    try (ObjectInputStream ois = new ObjectInputStream(GzipUtils.createGzipInputStream(new FileInputStream(pubFile)))) {
      key = (RSAPublicKey) ois.readObject();
    }
    try (FileOutputStream fos = new FileOutputStream(outFile)) {
      fromPublicKey(key, fos);
    }
    System.out.println("success");
  }
}
