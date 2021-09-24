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

import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import com.rtg.AbstractTest;
import com.rtg.util.StringUtils;
import com.rtg.util.io.TestDirectory;

/**
 * Test
 */
public class HashBlacklistTest extends AbstractTest {

  public void testBlacklistExists() throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      final File blacklistDir = new File(dir, HashBlacklist.BLACKLIST_SUBDIR);
      assertTrue(blacklistDir.mkdir());
      final File w15File = new File(blacklistDir, "w15");
      assertFalse(HashBlacklist.blacklistExists(dir, 15));
      assertTrue(w15File.createNewFile());
      assertTrue(HashBlacklist.blacklistExists(dir, 15));
    }
  }

  private static final String EXAMPLE_BLACKLIST = ""
    + "AAAAAAAA\t80" + StringUtils.LS
    + "AGAGAGAG\t100" + StringUtils.LS
    + "TTTTTTTT\t500" + StringUtils.LS
    + "CTCTCTCT\t1000" + StringUtils.LS;

  private static final Long[] EXPECTED_HASHES = {
    0b0000000000000000L,
    0b0010001000100010L,
    0b1111111111111111L,
    0b0111011101110111L
  };

  public void testLoad() throws IOException {
    List<Long> res = HashBlacklist.loadBlacklist(new StringReader(EXAMPLE_BLACKLIST), 0);
    assertEquals(Arrays.asList(EXPECTED_HASHES), res);
    res = HashBlacklist.loadBlacklist(new StringReader(EXAMPLE_BLACKLIST), 80);
    assertEquals(Arrays.asList(EXPECTED_HASHES), res);
    res = HashBlacklist.loadBlacklist(new StringReader(EXAMPLE_BLACKLIST), 81);
    assertEquals(Arrays.asList(EXPECTED_HASHES[1], EXPECTED_HASHES[2], EXPECTED_HASHES[3]), res);
    res = HashBlacklist.loadBlacklist(new StringReader(EXAMPLE_BLACKLIST), 101);
    assertEquals(Arrays.asList(EXPECTED_HASHES[2], EXPECTED_HASHES[3]), res);
    res = HashBlacklist.loadBlacklist(new StringReader(EXAMPLE_BLACKLIST), 501);
    assertEquals(Collections.singletonList(EXPECTED_HASHES[3]), res);
    res = HashBlacklist.loadBlacklist(new StringReader(EXAMPLE_BLACKLIST), 1001);
    assertEquals(Collections.emptyList(), res);
  }
}