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
package com.rtg.index;

import static com.rtg.util.StringUtils.LS;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import com.rtg.index.params.CreateParams;
import com.rtg.util.MapSet;
import com.rtg.util.PortableRandom;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.Talkback;
import com.rtg.util.test.NanoRegression;

import junit.framework.TestCase;

/**
 */
public abstract class AbstractIndexTest extends TestCase {

  protected NanoRegression mNano = null;

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
    mNano = new NanoRegression(this.getClass(), false);
  }

  @Override
  public void tearDown() throws Exception {
    // clear the module name so later tests don't report SlimException to the
    // Talkback system
    Talkback.setModuleName(null);
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
  }


  protected final boolean mCompressed;

  /**
   * @param compressed is compressed
   */
  public AbstractIndexTest(final boolean compressed) {
    mCompressed = compressed;
  }

  protected abstract IndexBase getIndex(final long size, final int hashBits, final Integer threshold);

  protected IndexBase getIndex(final long size, final int hashBits, final Integer threshold, final boolean twoPass) {
    if (twoPass) {
      return new IndexCompressed(new CreateParams(size, hashBits, hashBits, twoPass, false, false), threshold, false, threshold, threshold, 1);
    } else {
      return new IndexSimple(new CreateParams(size, hashBits, hashBits, twoPass, false, false), threshold, false, 1, threshold, threshold);
    }
  }

  protected void check(final IndexBase index, final long[] hashes, final int[] ids, final long[] missHashes, final boolean checkCount) throws IOException {
    check(null, index, hashes, ids, missHashes, checkCount);
  }

  protected abstract void check(String testId, final IndexBase index, final long[] hashes, final int[] ids, final long[] missHashes, final boolean checkCount) throws IOException;

  protected void checkEmpty(final IndexBase index, final long[] hashes, final boolean checkCount) throws IOException {
    final Finder finder = new Finder() {
      @Override
      public boolean found(final long id) {
        fail();
        return true;
      }
    };
    for (final long hash : hashes) {
      //System.err.println(" checkEmpty hash=" + hash);
      index.search(hash, finder);
      if (checkCount) {
        assertEquals(0, index.count(hash));
      }
      assertTrue(index.contains(hash) < 0);
    }
  }

  /**
   * Add the hashes and ids to the index.
   */
  protected void add(final IndexBase index, final long[] hashes, final int[] ids) {
    for (int i = 0; i < hashes.length; i++) {
      index.add(hashes[i], ids[i]);
    }
  }

  /** This tests the splitting and recombining of hashes into upper and lower bits */
  public void testCompressHash() {
    final IndexBase index = getIndex(1000000, 32, 20);
    for (long hash = 0; hash < 0x8800; hash += 0x1001) {
      final long upper = index.position(hash);
      final long lower = index.compressHash(hash);
      //      System.out.println(String.format("hash=%x, ", hash)
      //          + String.format("upper=%x, ", upper)
      //          + String.format("lower=%x", lower));
      assertEquals(hash, index.decompressHash(upper, lower));
    }
  }

  /**
   * Multiple entries for the same hash.
   * @throws IOException If an I/O error occurs
   */
  public void testOverflow2() throws IOException {
    final IndexBase hi = getIndex(10L, 16, 2);
    check("Overflow2", hi, new long[] {5, 1, 2, 5, 4, 2, 5, 5 }, new int[] {1, 2, 3, 4, 5, 6, 7, 8 }, new long[] {0, 3, 5, Short.MAX_VALUE - 1}, false);
  }

  /**
   * test to string empty
   * @throws IOException whenever.
   */
  public void testToStringEmpty() throws IOException {
    final IndexBase hi = getIndex(10L, 16, Integer.MAX_VALUE);
    mNano.check("ToStringEmpty1" + mCompressed, hi.toString() + hi.infoString());

    hi.freeze();
    if (mCompressed) {
      hi.freeze();
    }
    mNano.check("ToStringEmpty2" + mCompressed, hi.toString() + hi.perfString());
  }

  public void testSize1ToString() throws IOException {
    final IndexBase hi = getIndex(1, 16, 2);
    hi.add(0, 0);
    hi.freeze();
    if (mCompressed) {
      hi.add(0, 0);
      hi.freeze();
    }
    final String str = hi.toString();
    mNano.check("Size1ToString" + mCompressed, str);
  }

  /**
   * Empty index.
   */
  public void test0() throws IOException {
    final IndexBase hi = getIndex(1L, 16, Integer.MAX_VALUE);
    check(hi, new long[] {}, new int[] {}, new long[] {1, 0, 2, (1L << 16) - 1L}, true);
  }

  /**
   * Single entry.
   */
  public final void test1() throws IOException {
    final IndexBase hi = getIndex(10L, 16, Integer.MAX_VALUE);
    check(hi, new long[] {1 }, new int[] {1 }, new long[] {0, 2, Short.MAX_VALUE}, true);
  }

  /**
   * 16 bit hashes.
   */
  public final void testShort() throws IOException {
    final IndexBase hi = getIndex(10L, 16, Integer.MAX_VALUE);
    check(hi, new long[] {1, 2, 5, 42, Short.MAX_VALUE / 2, Short.MAX_VALUE }, new int[] {1, 2, 3, 4, 5, 6 }, new long[] {0, 3, Short.MAX_VALUE - 1}, true);
  }

  /**
   * 32 bt hashes.
   */
  public final void testInt() throws IOException {
    final IndexBase hi = getIndex(10L, 32, Integer.MAX_VALUE);
    check(hi, new long[] {1, 2, 5, 42, Integer.MAX_VALUE / 2, Integer.MAX_VALUE }, new int[] {1, 2, 3, 4, 5, 6 }, new long[] {0, 3, Short.MAX_VALUE - 1, Integer.MAX_VALUE - 1}, true);
  }

  /**
   * 64 bit hashes - this is important it winkled out problems with signed longs.
   */
  public final void testLong() throws IOException {
    final IndexBase hi = getIndex(10L, 64, Integer.MAX_VALUE);
    check("Long", hi, new long[] {-1, 1, 2, 5, 42, Long.MAX_VALUE / 2, Long.MAX_VALUE, Long.MIN_VALUE }, new int[] {1, 2, 3, 4, 5, 6, 7, 8 }, new long[] {Long.MIN_VALUE + 1, 0, 3, Short.MAX_VALUE - 1, Long.MAX_VALUE - 1}, true);
  }

  /**
   * Multiple entries for the same hash and duplicate ids for hash.
   */
  public final void testDuplicate() throws IOException {
    final IndexBase hi = getIndex(10L, 16, Integer.MAX_VALUE);
    check(hi, new long[] {5, 1, 2, 5, 4, 2, 5, 5 }, new int[] {1, 1, 1, 1, 1, 1, 1, 1 }, new long[] {0, 3, Short.MAX_VALUE - 1}, false);
  }

  /**
   * A 1000 random 16 bit entries.
   */
  public final void testBig1000() throws IOException {
    final int len = 1000;
    final long[] hashes = new long[len];
    final PortableRandom rand = new PortableRandom(42);
    final int[] ids = new int[len];
    for (int i = 0; i < len; i++) {
      hashes[i] = rand.nextInt(Short.MAX_VALUE + 1);
      ids[i] = rand.nextInt(Integer.MAX_VALUE);
    }
    final IndexBase hi = getIndex(len, 16, Integer.MAX_VALUE);
    check(hi, hashes, ids, new long[] {}, true);
  }

  /**
   * A 1000 random 16 bit entries chosen from a 20 member pool so that there are lots of duplicates.
   */
  public final void testBig1000pool() throws IOException {
    final int poolLen = 20;
    final long[] hpool = new long[poolLen];
    final int[] idpool = new int[poolLen];
    final PortableRandom rand = new PortableRandom(42);
    for (int i = 0; i < poolLen; i++) {
      hpool[i] = rand.nextInt(Short.MAX_VALUE + 1);
      idpool[i] = rand.nextInt(Integer.MAX_VALUE);
    }
    final int len = 1000;
    final long[] hashes = new long[len];
    final int[] ids = new int[len];
    for (int i = 0; i < len; i++) {
      hashes[i] = hpool[rand.nextInt(poolLen)];
      ids[i] = idpool[rand.nextInt(poolLen)];
    }

    final IndexBase hi = getIndex(len, 16, Integer.MAX_VALUE);
    check(hi, hashes, ids, new long[] {}, false);
  }

  public final void testDump() {
    checkDump(false);
    checkDump(true);
  }

  private void checkDump(final boolean twoPass) {
    final IndexBase hi = getIndex(8L, 16, Integer.MAX_VALUE, twoPass);
    for (int i = 0; i <= (twoPass ? 1 : 0); i++) {
      hi.add(1, 1);
      hi.add(1, 1);
      hi.add(2, 2);
      hi.add(2, 2);
      hi.add(2, 2);
      hi.add(3, 3);
      hi.freeze();
    }
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    final PrintStream pr = new PrintStream(bos);
    hi.dumpValues(pr);
    pr.close();
    final String exp = ""
        + "Index InitialPosition" + LS
        + "[0]  0" + LS
        + "[1]  6" + LS
        + "[2]  6" + LS
        + "[3]  6" + LS
        + "[4]  6" + LS
        + "[5]  6" + LS
        + "[6]  6" + LS
        + "[7]  6" + LS
        + "[8]  6" + LS
        + "[9]  6" + LS
        + "Index Hash  Values" + LS
        + "[0]  1  1" + LS
        + "[1]  1  1" + LS
        + "[2]  2  2" + LS
        + "[3]  2  2" + LS
        + "[4]  2  2" + LS
        + "[5]  3  3" + LS
        ;
    assertEquals(exp, bos.toString());
  }

  public final void testTrickySize1() {
    final IndexBase hi = getIndex(1, 16, 2);
    hi.freeze();
    if (mCompressed) {
      hi.freeze();
    }
    final String str = hi.toString();
    assertTrue(str.contains("pointer bits=2"));
  }

  public final void testTrickySize5() {
    final IndexBase hi = getIndex(5, 16, 2);
    hi.freeze();
    if (mCompressed) {
      hi.freeze();
    }
    assertTrue(hi.toString().contains("pointer bits=2"));
  }

  public final void testTrickySize9() {
    final IndexBase hi = getIndex(9, 16, 2);
    hi.freeze();
    if (mCompressed) {
      hi.freeze();
    }
    assertTrue(hi.toString().contains("pointer bits=3"));
  }

  public final void testThreshold0() throws IOException {
    checkThreshold0(false);
    checkThreshold0(true);
  }

  private void checkThreshold0(final boolean twoPass) throws IOException {
    final IndexBase hi = getIndex(10L, 16, 2, twoPass);
    for (int i = 0; i <= (twoPass ? 1 : 0); i++) {
      hi.add(1, 1);
      hi.add(1, 2);
      hi.add(1, 3);
      //System.err.println(hi.toString());
      hi.freeze();
    }
    //System.err.println(hi.toString());
    checkKeyCount(hi, 0, 0);
    checkKeyCount(hi, 1, 0);
    checkKeyCount(hi, 2, 0);
  }
  public void testMaxThreshold() throws IOException {
    final IndexBase hi = new IndexCompressed(new CreateParams(10L, 16, 16, true, false, false), 10, true, 3, 0, 1);
    for (int i = 0; i < 2; i++) {
      hi.add(1, 1);
      hi.add(1, 2);
      hi.add(1, 3);
      hi.add(2, 4);
      hi.add(2, 5);
      hi.add(0, 6);
      hi.add(0, 7);
      hi.add(0, 8);
      hi.add(0, 9);
      //System.err.println(hi.toString());
      hi.freeze();
    }
    //System.err.println(hi.toString());
    checkKeyCount(hi, 0, 0);
    checkKeyCount(hi, 1, 0);
    checkKeyCount(hi, 2, 2);
  }

  public final void testThreshold1() throws IOException {
    checkThreshold1(false);
    checkThreshold1(true);
  }

  public final void checkThreshold1(final boolean twoPass) throws IOException {
    final IndexBase hi = getIndex(10L, 16, 2, twoPass);
    for (int i = 0; i <= (twoPass ? 1 : 0); i++) {
      hi.add(0, 0);
      hi.add(1, 1);
      hi.add(1, 2);
      hi.add(1, 3);
      hi.add(2, 4);
      hi.add(2, 5);
      //System.err.println(hi.toString());
      hi.freeze();
    }
    //System.err.println(hi.toString());
    final Set<Integer> ss0 = searchSet(hi, 0);
    assertEquals(1, ss0.size());
    assertTrue(ss0.contains(0));

    final Set<Integer> ss1 = searchSet(hi, 1);
    assertEquals(0, ss1.size());

    final Set<Integer> ss2 = searchSet(hi, 2);
    assertEquals(2, ss2.size());
    assertTrue(ss2.contains(4));
    assertTrue(ss2.contains(5));

    final Set<Integer> ss3 = searchSet(hi, 3);
    assertEquals(0, ss3.size());

    checkKeyCount(hi, 0, 1);
    checkKeyCount(hi, 1, 0);
    checkKeyCount(hi, 2, 2);
  }

  public final void testThreshold2() throws IOException {
    checkThreshold2(false);
    checkThreshold2(true);
  }

  public final void checkThreshold2(final boolean twoPass) throws IOException {
    final IndexBase hi = getIndex(10L, 16, 2, twoPass);
    for (int i = 0; i <= (twoPass ? 1 : 0); i++) {
      hi.add(0, 0);
      hi.add(1, 1);
      hi.add(1, 2);
      hi.add(1, 3);
      hi.add(2, 4);
      hi.add(2, 5);
      hi.add(2, 6);
      //System.err.println(hi.toString());
      hi.freeze();
    }
    //System.err.println(hi.toString());
    final Set<Integer> ss0 = searchSet(hi, 0);
    assertEquals(1, ss0.size());
    assertTrue(ss0.contains(0));

    final Set<Integer> ss1 = searchSet(hi, 1);
    assertEquals(0, ss1.size());

    final Set<Integer> ss2 = searchSet(hi, 2);
    assertEquals(0, ss2.size());

    final Set<Integer> ss3 = searchSet(hi, 3);
    assertEquals(0, ss3.size());

    checkKeyCount(hi, 0, 1);
    checkKeyCount(hi, 1, 0);
    checkKeyCount(hi, 2, 0);
  }

  public final void testThreshold3() throws IOException {
    checkThreshold3(false);
    checkThreshold3(true);
  }

  public final void checkThreshold3(final boolean twoPass) throws IOException {
    final IndexBase hi = getIndex(10L, 16, 2, twoPass);
    for (int i = 0; i <= (twoPass ? 1 : 0); i++) {
      hi.add(0, 0);
      hi.add(1, 1);
      hi.add(1, 2);
      hi.add(1, 3);
      hi.add(2, 4);
      hi.add(2, 5);
      hi.add(2, 6);
      hi.add(3, 7);
      //System.err.println(hi.toString());
      hi.freeze();
    }
    //hi.dumpValues(System.out);
    final Set<Integer> ss0 = searchSet(hi, 0);
    assertEquals(1, ss0.size());
    assertTrue(ss0.contains(0));

    final Set<Integer> ss1 = searchSet(hi, 1);
    assertEquals(0, ss1.size());

    final Set<Integer> ss2 = searchSet(hi, 2);
    assertEquals(0, ss2.size());

    final Set<Integer> ss3 = searchSet(hi, 3);
    assertEquals(1, ss3.size());
    assertTrue(ss3.contains(7));
    //System.err.println(hi.toString());
    checkKeyCount(hi, 0, 1);
    checkKeyCount(hi, 1, 0);
    checkKeyCount(hi, 2, 0);
    checkKeyCount(hi, 3, 1);
  }

  public final void testThreshold4() throws IOException {
    checkThreshold4(false);
    checkThreshold4(true);
  }

  public final void checkThreshold4(final boolean twoPass) throws IOException {
    final IndexBase hi = getIndex(10L, 16, 2, twoPass);
    for (int i = 0; i <= (twoPass ? 1 : 0); i++) {
      hi.add(1, 1);
      hi.add(1, 2);
      hi.add(1, 3);
      hi.add(2, 4);
      hi.add(2, 5);
      hi.add(2, 6);
      hi.add(3, 7);
      //System.err.println(hi.toString());
      hi.freeze();
    }
    //System.err.println(hi.toString());


    final Set<Integer> ss0 = searchSet(hi, 0);
    assertEquals(0, ss0.size());

    final Set<Integer> ss1 = searchSet(hi, 1);
    assertEquals(0, ss1.size());

    final Set<Integer> ss2 = searchSet(hi, 2);
    assertEquals(0, ss2.size());

    final Set<Integer> ss3 = searchSet(hi, 3);
    assertEquals(1, ss3.size());
    assertTrue(ss3.contains(7));

    checkKeyCount(hi, 1, 0);
    checkKeyCount(hi, 2, 0);
    checkKeyCount(hi, 3, 1);
  }

  public final void testThresholdBug1() throws IOException {
    checkThresholdBug1(false);
    checkThresholdBug1(true);
  }

  public final void checkThresholdBug1(final boolean twoPass) throws IOException {
    final IndexBase hi = getIndex(10L, 16, 3, twoPass);
    for (int i = 0; i <= (twoPass ? 1 : 0); i++) {
      hi.add(1, 1);
      hi.add(1, 2);
      hi.add(1, 3);
      //System.err.println(hi.toString());
      hi.freeze();
    }
    //System.err.println(hi.toString());
    final Set<Integer> ss = searchSet(hi, 1);
    assertEquals(3, ss.size());
    assertTrue(ss.contains(1));
    assertTrue(ss.contains(2));
    assertTrue(ss.contains(3));

    checkKeyCount(hi, 1, 3);
  }

  private void checkKeyCount(final IndexBase hi, final long key, final int expected) throws IOException, IllegalStateException {
    final Set<Integer> actual = searchSet(hi, key);
    assertEquals(expected, actual.size());
    assertEquals(expected, hi.count(key));
  }

  private Set<Integer> searchSet(final IndexBase hi, final long key) throws IOException {
    final Set<Integer> actual = new HashSet<>();
    final Finder finder = new Finder() {
      @Override
      public boolean found(final long id) {
        actual.add((int) id);
        return true;
      }
    };
    hi.search(key, finder);
    return actual;
  }

  public final void testBadConstructor() {
    try {
      getIndex(-1L, 16, Integer.MAX_VALUE);
      fail();
    } catch (final RuntimeException e) {
      //expected
    }
    getIndex(1L, 1, Integer.MAX_VALUE);
    try {
      getIndex(1L, 0, Integer.MAX_VALUE);
      fail();
    } catch (final RuntimeException e) {
      //expected
    }
    try {
      getIndex(1L, 64, 0);
      fail();
    } catch (final RuntimeException e) {
      //expected
    }
  }

  protected void checkNoSearch(final IndexBase index) throws IOException {
    try {
      index.search(0, new Finder() {
        @Override
        public boolean found(final long id) {
          return true;
        }
      });
      fail();
    } catch (final IllegalStateException e) {
      //expected
    }
  }

  /*
   * Check that the hashes and ids can be found in the index.
   */
  protected void checkSearch(final IndexBase index, final long[] hashes, final int[] ids, final long[] missHashes, final boolean checkCount) throws IOException {
    final MapSet<Long, Long> ms = new MapSet<>();
    for (int i = 0; i < hashes.length; i++) {
      ms.put(hashes[i], (long) ids[i]);
    }
    for (final long missHashe : missHashes) {
      ms.remove(missHashe);
    }
    for (final Map.Entry<Long, Set<Long>> entry : ms.entrySet()) {
      final long hash = entry.getKey();
      final Set<Long> values = entry.getValue();
      final Set<Long> actual = new HashSet<>();
      final Finder finder = new Finder() {
        @Override
        public boolean found(final long id) {
          actual.add(id);
          return true;
        }
      };
      index.search(hash, finder);
      //assertEquals(values, actual);
      assertTrue(actual.equals(values));
      if (checkCount) {
        final int count = index.count(hash);
        //System.err.println("hash=" + hash + " count=" + count);
        assertEquals(actual.size(), count);
      }
    }
    for (final long hash : ms.keySet()) {
      final long found = index.contains(hash);
      assertTrue("hash=" + hash, found >= 0 && found < index.numberEntries());
      assertEquals(hash, index.getHash(found));
    }
    assertEquals(ms.numberOfKeys(), index.numberHashes());

    final MapSet<Long, Long> scanSet = new MapSet<>();
    index.scan(new FinderHashValue() {
      @Override
      public void found(long hash, long value) {
        scanSet.put(hash, value);
      }
    });
    assertEquals(ms, scanSet);

    //check search/getValue/getHash
    final MapSet<Long, Long> sv = new MapSet<>();
    for (long l = 0; l < index.numberEntries(); l++) {
      final long value = index.getValue(l);
      final long hash = index.getHash(l);
      sv.put(hash, value);
      final long find = index.contains(hash);
      assertEquals(hash, index.getHash(find));
    }
    assertEquals(ms, sv);
  }


}
