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

import java.net.NetworkInterface;
import java.net.SocketException;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.List;

/**
 * Determine the Ethernet address(es) of the current machine.
 *
 */
public class EthernetAddress {

  private final List<String> mAddresses;

  /**
   * @return the Ethernet address(es) of the current machine
   */
  public String[] getAddresses() {
    return mAddresses.toArray(new String[mAddresses.size()]);
  }

  /**
   * Determines the Ethernet address(es) of the current machine.
   *
   * @throws SocketException if an error occurred getting the Ethernet address(es)
   */
  public EthernetAddress() throws SocketException {
    mAddresses = new ArrayList<>();

    final Enumeration<NetworkInterface> interfaces = NetworkInterface.getNetworkInterfaces();
    while (interfaces.hasMoreElements()) {
      final NetworkInterface iface = interfaces.nextElement();
      final byte[] address = iface.getHardwareAddress();
      if (address == null || address.length != 6) {
        continue;
      }
      mAddresses.add(String.format("%02x %02x %02x %02x %02x %02x", address[0],
          address[1], address[2], address[3], address[4], address[5]));
    }
  }

  /**
   * Tests from the command line.
   *
   * @param args a <code>String</code> value
   * @throws SocketException if an error occurred getting the Ethernet address(es)
   */
  public static void main(final String[] args) throws SocketException {
    final EthernetAddress ea = new EthernetAddress();

    System.out.println("OS: " + System.getProperties().getProperty("os.name"));
    System.out.println("CPU cores: " + Runtime.getRuntime().availableProcessors());
    System.out.println("Mac Addresses: ");
    for (String s : ea.getAddresses()) {
      System.out.println("  " + s);
    }
  }
}

