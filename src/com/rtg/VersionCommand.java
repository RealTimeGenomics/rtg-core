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
package com.rtg;

import static com.rtg.util.StringUtils.LS;

import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintStream;

import com.rtg.util.Constants;
import com.rtg.util.Environment;
import com.rtg.util.License;
import com.rtg.util.diagnostic.SlimException;
import com.rtg.util.diagnostic.Talkback;

/**
 * Provides version information on <code>SLIM</code>
 */
public final class VersionCommand {

  private VersionCommand() { }


  static final String PROGRAM_STR = "Product:";
  static final String VERSION_STR = "Core Version:";
  static final String CONTACT_STR = "Contact:";

  static final String RAM_STR = "RAM:";
  static final String JAVA_VERSION_STR = "JVM:";

  static final String COPYRIGHT_STR = "(c) Real Time Genomics, 2014";

  static final String CITE_STR = "Citation:" + LS
    + "John G. Cleary, Ross Braithwaite, Kurt Gaastra, Brian S. Hilbush, Stuart Inglis, "
    + "Sean A. Irvine, Alan Jackson, Richard Littin, Sahar Nohzadeh-Malakshah, "
    + "Mehul Rathod, David Ware, Len Trigg, and Francisco M. De La Vega. "
    + "\"Joint Variant and De Novo Mutation Identification on Pedigrees from High-Throughput Sequencing Data.\" "
    + "Journal of Computational Biology. June 2014, 21(6): 405-419. doi:10.1089/cmb.2014.0029.";

  static final String PATENTS_STR = "Patents / Patents pending:" + LS
    + "US: 7,640,256, 13/129,329, 13/681,046, 13/681,215, 13/848,653, 13/925,704, 14/015,295, 13/971,654, 13/971,630, 14/564,810" + LS
    + "UK: 1222923.3, 1222921.7, 1304502.6, 1311209.9, 1314888.7, 1314908.3" + LS
    + "New Zealand: 626777, 626783, 615491, 614897, 614560" + LS
    + "Australia: 2005255348, Singapore: 128254";

  /**
   * Main function, entry-point for help
   *
   * @param args additional command-line arguments
   * @param out regular output stream, only for passing to delegate help printers
   * @return shell return code 0 for success, anything else for failure
   */
  public static int mainInit(String[] args, final OutputStream out) {
    final OutputStreamWriter wr = new OutputStreamWriter(out);
    try {
      try {
        wr.write(getVersion());

        if (args.length > 0 && "--environment".equals(args[0])) {
          wr.write(Talkback.getEnvironment());
        }
      } finally {
        wr.flush();
      }
    } catch (final IOException e) {
      throw new SlimException(e);
    }
    return 0;
  }

  /**
   * Main function, entry-point for help
   * These next two methods are strictly for sharpen, do not remove them however stupid they look.
   * @param out regular output stream, only for passing to delegate help printers
   * @return shell return code 0 for success, anything else for failure
   */
  public static int mainInit(final PrintStream out) {
    printVersion(out);
    return 0;
  }
  private static void printVersion(final PrintStream out) {
    try {
      out.print(getVersion());
    } finally {
      out.flush();
    }
  }


  static String getRamString() {
    final double maxMemory = Runtime.getRuntime().maxMemory() * 10.0 / 1024 / 1024 / 1024;
    final double totalMemory = Environment.getTotalMemory() * 10.0 / 1024 / 1024 / 1024;
    final int percentageAllocated = (int) (maxMemory * 100.0 / totalMemory);
    return getRamString(percentageAllocated, maxMemory, totalMemory);
  }

  static String getRamString(final int percentageAllocated, final double maxMemory, final double totalMemory) {
    return ((int) maxMemory) / 10.0 + "GB of " + ((int) totalMemory) / 10.0 + "GB RAM can be used by rtg "
    + "(" + percentageAllocated + "%)";
  }

  static String getJavaVersion() {
    final String version = System.getProperty("java.version");
    final String vm = System.getProperty("java.vm.name");
    final StringBuilder sb = new StringBuilder();
    if (vm != null) {
      sb.append(vm);
    }
    if (version != null) {
      if (sb.length() != 0) {
        sb.append(" ");
      }
      sb.append(version);
    }
    if (sb.length() == 0) {
      return "unknown";
    }
    return sb.toString();
  }

  /**
   * Returns the version string.
   * @return the version string.
   */
  public static String getVersion() {
    final StringBuilder sb = new StringBuilder();
    sb.append(PROGRAM_STR + " ").append(Environment.getProductName()).append(LS);
    sb.append(VERSION_STR + " ").append(Environment.getCoreVersion()).append(LS);
    if (License.checkLicense()) {
      sb.append(RAM_STR).append(" ").append(getRamString()).append(LS);
      sb.append(JAVA_VERSION_STR).append(" ").append(getJavaVersion()).append(LS);
    }
    sb.append(LicenseCommand.getLicenseSummary());
    sb.append(CONTACT_STR + " " + Constants.SUPPORT_EMAIL_ADDR).append(LS);
    sb.append(LS);
    sb.append(PATENTS_STR).append(LS);
    sb.append(LS);
    sb.append(CITE_STR).append(LS);
    sb.append(LS);
    sb.append(COPYRIGHT_STR).append(LS);
    sb.append(LS);
    return sb.toString();
  }

}
