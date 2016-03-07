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
package com.rtg.index.params;

import com.rtg.index.HashBitHandle;
import com.rtg.index.IndexBase;
import com.rtg.index.IndexSimple;
import com.rtg.index.IndexUtils;
import com.rtg.index.RepeatFrequencyFilterMethod;
import com.rtg.index.UnfilteredFilterMethod;
import com.rtg.index.params.CreateParams.CreateParamsBuilder;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.array.ArrayHandle;
import com.rtg.util.array.ArrayType;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 */
public class CreateParamsTest extends TestCase {

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }

  CreateParams getParams(final long size, final int hashBits) {
    return getParams(size, hashBits, hashBits);
  }

  CreateParams getParams(final long size, final int hashBits, final int windowBits) {
    return new CreateParams(size, hashBits, windowBits, 31, false, true, false, false);
  }

  public void testEquals() throws Exception {
    final CreateParams a1 = getParams(100, 2);
    final CreateParams a2 = getParams(100, 2);
    final CreateParams b = getParams(101, 2);
    final CreateParams c = getParams(100, 3);
    TestUtils.equalsHashTest(new CreateParams[][] {{a1, a2}, {b}, {c}});
  }

  public void testWindows1() {
    final CreateParams ip = getParams(0L, 64, 64);
    ip.integrity();
    assertEquals(0, ip.size());
    assertEquals(64, ip.hashBits());
    assertEquals(64, ip.windowBits());
    final String str = ip.toString();
    assertEquals(" size=0 hash bits=64 initial pointer bits=2 value bits=31", str);
  }

  public void testWindows2() {
    final CreateParams ip = getParams(0L, 64, 65);
    ip.integrity();
    assertEquals(0, ip.size());
    assertEquals(64, ip.hashBits());
    assertEquals(65, ip.windowBits());
    final String str = ip.toString();
    assertEquals(" size=0 hash bits=64 initial pointer bits=2 value bits=31", str);
  }

  private static final String EXPECTED0 = ""
      + StringUtils.LS
      + "Memory Usage\tbytes\tlength" + StringUtils.LS
      + "\t\t1\t1\tHash" + StringUtils.LS
      + "\t\t4\t1\tValue" + StringUtils.LS
      + "\t\t6\t6\tInitial Position" + StringUtils.LS
      + "\t\t4\t16\tBit vector" + StringUtils.LS
      + "\t\t15\t\tTotal" + StringUtils.LS
      ;

  private static final String MEM_EXPECTED = ""
      + "\tMemory\tHash\t1" + StringUtils.LS
      + "\tMemory\tValue\t4" + StringUtils.LS
      + "\tMemory\tInitial_position\t6" + StringUtils.LS
      + "\tMemory\tBit_vector\t4" + StringUtils.LS;

  private static final String HASH_BUCKET = ""
      + "Hash counts\t0\t1\t2" + StringUtils.LS
      + "\t\t0\t0\t0" + StringUtils.LS
      + "Bucket counts\t0\t1\t2" + StringUtils.LS
      + "\t\t0\t0\t0" + StringUtils.LS
      ;

  public void test0() throws Exception {
    final CreateParams ip = getParams(1L, 8);
    ip.integrity();
    assertEquals(1, ip.size());
    assertEquals(8, ip.hashBits());
    assertEquals(8, ip.windowBits());
    final String str = ip.toString();
    assertEquals(" size=1 hash bits=8 initial pointer bits=2 value bits=31", str);

    //test all the handle methods
    final ArrayHandle ha = ip.hash();
    assertEquals("BYTE", ha.type().toString());
    assertEquals(1, ha.length());
    assertEquals(1, ha.bytes());

    final ArrayHandle ia = ip.value();
    assertEquals(ArrayType.INTEGER , ia.type());
    assertEquals(1, ia.length());
    assertEquals(4, ia.bytes());

    assertEquals(2, ip.initialPointerBits());
    final ArrayHandle pa = ip.initialPosition();
    assertEquals("BYTE" , pa.type().toString());
    assertEquals(6, pa.length());
    assertEquals(6, pa.bytes());

    final HashBitHandle hbh = ip.bitVector();
    assertEquals(16, hbh.length());
    assertEquals(4, hbh.bytes());
    assertEquals(4, hbh.bits());

    assertEquals(15, IndexUtils.bytes(ip));
    assertEquals(EXPECTED0, IndexUtils.memString(ip));
    assertEquals(MEM_EXPECTED, IndexUtils.memToString(ip));


    final IndexBase ii = new IndexSimple(ip, new UnfilteredFilterMethod(), 1);
    assertEquals(15, ii.bytes());
    assertEquals(EXPECTED0 + HASH_BUCKET, ii.infoString());
  }

  private static final String EXPECTED42 = ""
      + StringUtils.LS
      + "Memory Usage\tbytes\tlength" + StringUtils.LS
      + "\t\t42\t42\tHash" + StringUtils.LS
      + "\t\t168\t42\tValue" + StringUtils.LS
      + "\t\t34\t34\tInitial Position" + StringUtils.LS
      + "\t\t32\t256\tBit vector" + StringUtils.LS
      + "\t\t276\t\tTotal" + StringUtils.LS;

  public void test42() throws Exception {
    final CreateParams ip = getParams(42L, 8);
    ip.integrity();
    assertEquals(42, ip.size());
    assertEquals(8, ip.hashBits());
    assertEquals(8, ip.windowBits());
    final String str = ip.toString();
    assertEquals(" size=42 hash bits=8 initial pointer bits=5 value bits=31", str);

    //test all the handle methods
    final ArrayHandle ha = ip.hash();
    assertEquals("BYTE", ha.type().toString());
    assertEquals(42, ha.length());
    assertEquals(42, ha.bytes());

    final ArrayHandle ia = ip.value();
    assertEquals(ArrayType.INTEGER , ia.type());
    assertEquals(42, ia.length());
    assertEquals(168, ia.bytes());

    assertEquals(5, ip.initialPointerBits());
    final ArrayHandle pa = ip.initialPosition();
    assertEquals("BYTE" , pa.type().toString());
    assertEquals(34, pa.length());
    assertEquals(34, pa.bytes());

    final HashBitHandle hbh = ip.bitVector();
    assertEquals(256, hbh.length());
    assertEquals(32, hbh.bytes());
    assertEquals(8, hbh.bits());

    assertEquals(276, IndexUtils.bytes(ip));
    assertEquals(EXPECTED42, IndexUtils.memString(ip));

    final IndexBase ii = new IndexSimple(ip, new UnfilteredFilterMethod(), 1);
    assertEquals(276, ii.bytes());
    assertEquals(EXPECTED42 + HASH_BUCKET, ii.infoString());
  }

  private static final String EXPECTED42_A = ""
      + StringUtils.LS
      + "Memory Usage\tbytes\tlength" + StringUtils.LS
      + "\t\t42\t42\tHash" + StringUtils.LS
      + "\t\t168\t42\tValue" + StringUtils.LS
      + "\t\t34\t34\tInitial Position" + StringUtils.LS
      + "\t\t8\t64\tBit vector" + StringUtils.LS
      + "\t\t252\t\tTotal" + StringUtils.LS;

  public void test42a() throws Exception {
    final CreateParams ip = getParams(42L, 1);
    ip.integrity();
    assertEquals(42, ip.size());
    assertEquals(1, ip.hashBits());
    assertEquals(1, ip.windowBits());
    final String str = ip.toString();
    assertEquals(" size=42 hash bits=1 initial pointer bits=5 value bits=31", str);

    final HashBitHandle hbh = ip.bitVector();
    assertEquals(64, hbh.length());
    assertEquals(8, hbh.bytes());
    assertEquals(6, hbh.bits());

    assertEquals(252, IndexUtils.bytes(ip));
    assertEquals(EXPECTED42_A, IndexUtils.memString(ip));

    final IndexBase ii = new IndexSimple(ip, new RepeatFrequencyFilterMethod(1000, false, Integer.MAX_VALUE, 0), 1);
    assertEquals(252, ii.bytes());
    assertEquals(EXPECTED42_A + HASH_BUCKET, ii.infoString());
  }

  private static final String EXPECTED42_B = ""
      + StringUtils.LS
      + "Memory Usage\tbytes\tlength" + StringUtils.LS
      + "\t\t42\t42\tHash" + StringUtils.LS
      + "\t\t168\t42\tValue" + StringUtils.LS
      + "\t\t34\t34\tInitial Position" + StringUtils.LS
      + "\t\t8\t64\tBit vector" + StringUtils.LS
      + "\t\t252\t\tTotal" + StringUtils.LS;

  public void test42b() throws Exception {
    final CreateParams ip = getParams(42L, 2);
    ip.integrity();
    assertEquals(42, ip.size());
    assertEquals(2, ip.hashBits());
    assertEquals(2, ip.windowBits());
    final String str = ip.toString();
    assertEquals(" size=42 hash bits=2 initial pointer bits=5 value bits=31", str);

    final HashBitHandle hbh = ip.bitVector();
    assertEquals(64, hbh.length());
    assertEquals(8, hbh.bytes());
    assertEquals(6, hbh.bits());

    assertEquals(252, IndexUtils.bytes(ip));
    assertEquals(EXPECTED42_B, IndexUtils.memString(ip));

    final IndexBase ii = new IndexSimple(ip, new RepeatFrequencyFilterMethod(1000, false, Integer.MAX_VALUE, 0), 1);
    assertEquals(252, ii.bytes());
    assertEquals(EXPECTED42_B + HASH_BUCKET, ii.infoString());
  }

  private static final String EXPECTED42_C = ""
      + StringUtils.LS
      + "Memory Usage\tbytes\tlength" + StringUtils.LS
      + "\t\t42\t42\tHash" + StringUtils.LS
      + "\t\t168\t42\tValue" + StringUtils.LS
      + "\t\t34\t34\tInitial Position" + StringUtils.LS
      + "\t\t8\t64\tBit vector" + StringUtils.LS
      + "\t\t252\t\tTotal" + StringUtils.LS;

  public void test42c() throws Exception {
    final CreateParams ip = getParams(42L, 3);
    ip.integrity();
    assertEquals(42, ip.size());
    assertEquals(3, ip.hashBits());
    assertEquals(3, ip.windowBits());
    final String str = ip.toString();
    assertEquals(" size=42 hash bits=3 initial pointer bits=5 value bits=31", str);

    final HashBitHandle hbh = ip.bitVector();
    assertEquals(64, hbh.length());
    assertEquals(8, hbh.bytes());
    assertEquals(6, hbh.bits());

    assertEquals(252, IndexUtils.bytes(ip));
    assertEquals(EXPECTED42_C, IndexUtils.memString(ip));

    final IndexBase ii = new IndexSimple(ip, new RepeatFrequencyFilterMethod(1000, false, Integer.MAX_VALUE, 0), 1);
    assertEquals(252, ii.bytes());
    assertEquals(EXPECTED42_C + HASH_BUCKET, ii.infoString());
  }

  private static final String EXPECTED42_D = ""
      + StringUtils.LS
      + "Memory Usage\tbytes\tlength" + StringUtils.LS
      + "\t\t42\t42\tHash" + StringUtils.LS
      + "\t\t168\t42\tValue" + StringUtils.LS
      + "\t\t34\t34\tInitial Position" + StringUtils.LS
      + "\t\t8\t64\tBit vector" + StringUtils.LS
      + "\t\t252\t\tTotal" + StringUtils.LS;

  public void test42d() throws Exception {
    final CreateParams ip = getParams(42L, 4);
    ip.integrity();
    assertEquals(42, ip.size());
    assertEquals(4, ip.hashBits());
    assertEquals(4, ip.windowBits());
    final String str = ip.toString();
    assertEquals(" size=42 hash bits=4 initial pointer bits=5 value bits=31", str);

    final HashBitHandle hbh = ip.bitVector();
    assertEquals(64, hbh.length());
    assertEquals(8, hbh.bytes());
    assertEquals(6, hbh.bits());

    assertEquals(252, IndexUtils.bytes(ip));
    assertEquals(EXPECTED42_D, IndexUtils.memString(ip));

    final IndexBase ii = new IndexSimple(ip, new RepeatFrequencyFilterMethod(1000, false, Integer.MAX_VALUE, 0), 1);
    assertEquals(252, ii.bytes());
    assertEquals(EXPECTED42_D + HASH_BUCKET, ii.infoString());
  }

  private static final String EXPECTED42_E = ""
      + StringUtils.LS
      + "Memory Usage\tbytes\tlength" + StringUtils.LS
      + "\t\t42\t42\tHash" + StringUtils.LS
      + "\t\t168\t42\tValue" + StringUtils.LS
      + "\t\t34\t34\tInitial Position" + StringUtils.LS
      + "\t\t16\t128\tBit vector" + StringUtils.LS
      + "\t\t260\t\tTotal" + StringUtils.LS;

  public void test42e() throws Exception {
    final CreateParams ip = getParams(42L, 7);
    ip.integrity();
    assertEquals(42, ip.size());
    assertEquals(7, ip.hashBits());
    assertEquals(7, ip.windowBits());
    final String str = ip.toString();
    assertEquals(" size=42 hash bits=7 initial pointer bits=5 value bits=31", str);

    final HashBitHandle hbh = ip.bitVector();
    assertEquals(7, hbh.bits());
    assertEquals(128, hbh.length());
    assertEquals(16, hbh.bytes());

    assertEquals(260, IndexUtils.bytes(ip));
    assertEquals(EXPECTED42_E, IndexUtils.memString(ip));

    final IndexBase ii = new IndexSimple(ip, new RepeatFrequencyFilterMethod(1000, false, Integer.MAX_VALUE, 0), 1);
    assertEquals(260, ii.bytes());
    assertEquals(EXPECTED42_E + HASH_BUCKET, ii.infoString());
  }

  private static final String EXPECTED42_F = ""
      + StringUtils.LS
      + "Memory Usage\tbytes\tlength" + StringUtils.LS
      + "\t\t120\t15\tHash" + StringUtils.LS
      + "\t\t60\t15\tValue" + StringUtils.LS
      + "\t\t10\t10\tInitial Position" + StringUtils.LS
      + "\t\t32\t256\tBit vector" + StringUtils.LS
      + "\t\t222\t\tTotal" + StringUtils.LS;

  public void test42f() throws Exception {
    final CreateParams ip = getParams(15L, 64);
    ip.integrity();
    assertEquals(15, ip.size());
    assertEquals(64, ip.hashBits());
    assertEquals(64, ip.windowBits());
    final String str = ip.toString();
    assertEquals(" size=15 hash bits=64 initial pointer bits=3 value bits=31", str);

    final HashBitHandle hbh = ip.bitVector();
    assertEquals(8, hbh.bits());
    assertEquals(256, hbh.length());
    assertEquals(32, hbh.bytes());

    assertEquals(222, IndexUtils.bytes(ip));
    assertEquals(EXPECTED42_F, IndexUtils.memString(ip));

    final IndexBase ii = new IndexSimple(ip, new RepeatFrequencyFilterMethod(1000, false, Integer.MAX_VALUE, 0), 1);
    assertEquals(222, ii.bytes());
    assertEquals(EXPECTED42_F  + HASH_BUCKET, ii.infoString());
  }

  private static final String EXPECTED42_G = ""
      + StringUtils.LS
      + "Memory Usage\tbytes\tlength" + StringUtils.LS
      + "\t\t64\t8\tHash" + StringUtils.LS
      + "\t\t32\t8\tValue" + StringUtils.LS
      + "\t\t10\t10\tInitial Position" + StringUtils.LS
      + "\t\t16\t128\tBit vector" + StringUtils.LS
      + "\t\t122\t\tTotal" + StringUtils.LS;

  public void test42g() throws Exception {
    final CreateParams ip = getParams(8L, 64);
    ip.integrity();
    assertEquals(8, ip.size());
    assertEquals(64, ip.hashBits());
    assertEquals(64, ip.windowBits());
    final String str = ip.toString();
    assertEquals(" size=8 hash bits=64 initial pointer bits=3 value bits=31", str);

    final HashBitHandle hbh = ip.bitVector();
    assertEquals(7, hbh.bits());
    assertEquals(128, hbh.length());
    assertEquals(16, hbh.bytes());

    assertEquals(122, IndexUtils.bytes(ip));
    assertEquals(EXPECTED42_G, IndexUtils.memString(ip));

    final IndexBase ii = new IndexSimple(ip, new RepeatFrequencyFilterMethod(1000, false, Integer.MAX_VALUE, 0), 1);
    assertEquals(122, ii.bytes());
    assertEquals(EXPECTED42_G + HASH_BUCKET, ii.infoString());
  }

  private static final String EXPECTED42_Z = ""
      + StringUtils.LS
      + "Memory Usage\tbytes\tlength" + StringUtils.LS
      + "\t\t336\t42\tHash" + StringUtils.LS
      + "\t\t168\t42\tValue" + StringUtils.LS
      + "\t\t34\t34\tInitial Position" + StringUtils.LS
      + "\t\t128\t1,024\tBit vector" + StringUtils.LS
      + "\t\t666\t\tTotal" + StringUtils.LS;

  public void test42z() throws Exception {
    final CreateParams ip = getParams(42L, 64);
    ip.integrity();
    assertEquals(42, ip.size());
    assertEquals(64, ip.hashBits());
    final String str = ip.toString();
    assertEquals(" size=42 hash bits=64 initial pointer bits=5 value bits=31", str);

    final HashBitHandle hbh = ip.bitVector();
    assertEquals(10, hbh.bits());
    assertEquals(1024, hbh.length());
    assertEquals(128, hbh.bytes());

    assertEquals(666, IndexUtils.bytes(ip));
    assertEquals(EXPECTED42_Z, IndexUtils.memString(ip));

    final IndexBase ii = new IndexSimple(ip, new RepeatFrequencyFilterMethod(1000, false, Integer.MAX_VALUE, 0), 1);
    assertEquals(666, ii.bytes());
    assertEquals(EXPECTED42_Z + HASH_BUCKET, ii.infoString());
  }

  private static final long SIZE = (1L << 29) + 2;
  private static final long PLENGTH = (1L << 29) + 2;

  public void testBig() throws Exception {
    final CreateParams ip = getParams(SIZE, 8);
    ip.integrity();
    assertEquals(SIZE, ip.size());
    assertEquals(8, ip.hashBits());
    assertEquals(8, ip.windowBits());
    final String str = ip.toString();
    assertEquals(" size=536,870,914 hash bits=8 initial pointer bits=29 value bits=31", str);

    //test all the handle methods
    final ArrayHandle ha = ip.hash();
    assertEquals("BYTE", ha.type().toString());
    assertEquals(SIZE, ha.length());
    assertEquals(536870914, ha.bytes());

    final ArrayHandle ia = ip.value();
    assertEquals(ArrayType.INTEGER , ia.type());
    assertEquals(SIZE, ia.length());
    assertEquals(4 * SIZE, ia.bytes());

    assertEquals(29, ip.initialPointerBits());
    final ArrayHandle pa = ip.initialPosition();
    assertEquals(ArrayType.INTEGER , pa.type());
    assertEquals(PLENGTH, pa.length());
    assertEquals(4 * PLENGTH, pa.bytes());

    final HashBitHandle hbh = ip.bitVector();
    assertEquals(1073741824, hbh.length());
    assertEquals(134217728, hbh.bytes());
    assertEquals(30, hbh.bits());

    assertEquals(4966055954L, IndexUtils.bytes(ip));
  }

  public void testComputeInitialPointerBitsIdeal1() {
    assertEquals(2, CreateParams.computeInitialPointerBitsIdeal(0, 0, false));
    assertEquals(24, CreateParams.computeInitialPointerBitsIdeal(4999999981L, 40, true));
    assertEquals(24, CreateParams.computeInitialPointerBitsIdeal(4999999981L, 40, false));
    assertEquals(28, CreateParams.computeInitialPointerBitsIdeal(4999999981L, 36, true));
    assertEquals(10, CreateParams.computeInitialPointerBitsIdeal(68830001L, 64, true));
    final ArrayHandle handle = getInitialPositionHandle(68830001L, 64);
    assertEquals(4104, handle.bytes());
    assertEquals(1026, handle.length());
  }

  public void testComputeInitialPointerBitsIdeal2() {
    for (int h = 0; h <= 64; h += 5) {
      int pl = CreateParams.computeInitialPointerBitsIdeal(1, h, true);
      for (long l = 1; l <= Long.MAX_VALUE >> 2; l += l) {
        final int p = CreateParams.computeInitialPointerBitsIdeal(l, h, true);
        //System.err.println("c=" + c + " l=" + l + " p=" + p);
        assertTrue(p >= pl);
        pl = p;
      }
    }
  }

  private ArrayHandle getInitialPositionHandle(long size, int hashbits) {
    return new CreateParams(size, hashbits, hashbits, 31, false, true, true, true).initialPosition();
  }

  public void testBuilder() {

    final CreateParams def = new CreateParamsBuilder().create();
    assertEquals(20, def.size());
    assertEquals(20, def.hashBits());
    assertEquals(20, def.windowBits());
    assertEquals(20, def.valueBits());
    assertTrue(def.compressHashes());
    assertNotNull(def.bitVector());

    final CreateParamsBuilder b = new CreateParamsBuilder();
    assertEquals(b, b.compressHashes(true));
    assertEquals(b, b.createBitVector(false));
    assertEquals(b, b.spaceEfficientButUnsafe(true));
    assertEquals(b, b.hashBits(35));
    assertEquals(b, b.size(18703));
    assertEquals(b, b.valueBits(57));
    assertEquals(b, b.windowBits(35));
    final CreateParams p = b.create();

    assertTrue(p.compressHashes());
    assertNull(p.bitVector());
    assertEquals(35, p.hashBits());
    assertEquals(18703, p.size());
    assertEquals(57, p.valueBits());
    assertEquals(35, p.windowBits());
    assertEquals(14, p.initialPointerBits());

    assertEquals(b, b.ideal(true));
    final CreateParams p2 = b.create();
    assertEquals(10, p2.initialPointerBits());
  }

  //test that will work ok when using the extended interface in IndexCompressed
  public void testExtendedLoop() {
    for (int i = 1; i <= 94; i++) {
      final CreateParams params = new CreateParams(Integer.MAX_VALUE, i, i, 0, true, true, true, false);
      params.integrity();
      assertEquals(Math.max(0, i - 30), params.hashCompressedBits());
    }
  }

  //test that will work ok when using the extended interface in IndexCompressed
  public void testExtendedLoop2() {
    for (int i = 64; i <= 94; i++) {
      final CreateParams params = new CreateParams(Integer.MAX_VALUE, i, i, 0, true, false, true, true);
      params.integrity();
      assertEquals(Math.min(64, i - 10), params.hashCompressedBits());
      //System.err.println(i + ":" + params.hashCompressedBits());
    }
  }


  public void testHashBits() {
    assertEquals(1, CreateParams.calculateHashBits(1, 1));
    assertEquals(2, CreateParams.calculateHashBits(2, 1));
    assertEquals(2, CreateParams.calculateHashBits(1, 2));
    assertEquals(62, CreateParams.calculateHashBits(2, 31));
    assertEquals(64, CreateParams.calculateHashBits(2, 32));
    assertEquals(64, CreateParams.calculateHashBits(2, 33));
    assertEquals(60, CreateParams.calculateHashBits(5, 12));
    assertEquals(64, CreateParams.calculateHashBits(5, 13));
  }
}

