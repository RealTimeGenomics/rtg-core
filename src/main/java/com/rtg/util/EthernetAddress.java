/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package com.rtg.util;

import java.net.NetworkInterface;
import java.net.SocketException;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.List;

/**
 * Determine the Ethernet address(es) of the current machine.
 */
public class EthernetAddress {

  private final List<String> mAddresses;

  /**
   * @return the Ethernet address(es) of the current machine
   */
  public String[] getAddresses() {
    return mAddresses.toArray(new String[0]);
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

    System.out.println("OS: " + System.getProperty("os.name"));
    System.out.println("CPU cores: " + Runtime.getRuntime().availableProcessors());
    System.out.println("Mac Addresses: ");
    for (String s : ea.getAddresses()) {
      System.out.println("  " + s);
    }
  }
}

