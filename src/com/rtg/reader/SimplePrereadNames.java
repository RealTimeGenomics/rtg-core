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
package com.rtg.reader;

import java.io.IOException;
import java.io.OutputStream;

import com.rtg.util.array.byteindex.ByteChunks;
import com.rtg.util.array.longindex.LongChunks;

/**
 * Simple implementation of {@link PrereadNames}
 */
public class SimplePrereadNames implements PrereadNamesInterface {

  private final ByteChunks mNameBytes = new ByteChunks(0);
  private final LongChunks mPointers = new LongChunks(1);
  private long mTotalNameSize = 0;

  /**
   * Default constructor
   */
  public SimplePrereadNames() {
    //This is here to make the first pointer explicit
    mPointers.set(0, 0);
  }

  @Override
  public long length() {
    return mPointers.length() - 1;
  }

  @Override
  public String name(long id) {
    return new String(getNameBytes(id));
  }

  byte[] getNameBytes(long id) {
    final long start = mPointers.get(id);
    final long end = mPointers.get(id + 1);
    final int len = (int) (end - start);
    final byte[] out = new byte[len];
    mNameBytes.getBytes(out, 0, start, len);
    return out;
  }

  /**
   * Add a name after the given id
   * @param id please make this the current length of this
   * @param name the name to add
   */
  public void setName(long id, String name) {
    assert id == length();
    final byte[] nameBytes = name.getBytes();
    final int length = nameBytes.length;
    mNameBytes.extendBy(length);
    mNameBytes.copyBytes(nameBytes, 0, mTotalNameSize, length);
    mTotalNameSize += length;
    mPointers.extendBy(1);
    mPointers.set(id + 1, mTotalNameSize);
  }

  /**
   * Calculate the checksum of the names in a manner compatible with
   * how the checksum is calculated in the SDF.
   *
   * @return the checksum of the names.
   */
  @Override
  public long calcChecksum() {
    final PrereadHashFunction namef = new PrereadHashFunction();
    for (long k = 0; k < length(); k++) {
      final byte[] name = getNameBytes(k);
      namef.irvineHash(name);
      namef.irvineHash(name.length);
    }
    return namef.getHash();
  }

  /**
   * Returns size of object in bytes
   * @return size of object in no of bytes
   */
  @Override
  public long bytes() {
    //low balling estimate at 2 pointers and a long per entry. TODO make this more reasonable
    return mNameBytes.bytes() + mPointers.bytes();
  }

  @Override
  public void writeName(Appendable a, long id) throws IOException {
    a.append(name(id));
  }

  @Override
  public void writeName(OutputStream stream, long id) throws IOException {
    stream.write(getNameBytes(id));
  }

}
