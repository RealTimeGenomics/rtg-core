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

import java.net.SocketException;
import java.util.regex.Pattern;

import junit.framework.TestCase;

/**
 */
public class EthernetAddressTest extends TestCase {

  public void testEthernetAddress() throws SocketException {
    final EthernetAddress ea = new EthernetAddress();
    final String[] addresses = ea.getAddresses();
    assertNotNull(addresses);
    assertTrue(addresses.length > 0);

    for (final String a : addresses) {
      assertTrue(Pattern.compile("([0-9a-f]{2} ){5}[0-9a-f]{2}").matcher(a).matches());
    }
  }
}

