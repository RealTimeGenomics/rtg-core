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
package com.rtg.usage;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Properties;
import java.util.UUID;

import com.rtg.util.Environment;
import com.rtg.util.License;
import com.rtg.util.MD5Utils;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;

import junit.framework.TestCase;

/**
 */
public class UsageServerCliTest extends TestCase {

  public void test() throws Exception {
    try (TestDirectory dir = new TestDirectory(); TestDirectory configDir = new TestDirectory()) {
      final File configFile = UsageConfiguration.createSimpleConfigurationFile(new File(configDir, "config"), dir.getPath(), null);
      final UsageServerCli cli = new UsageServerCli() {
        @Override
        protected int mainExec(OutputStream out, PrintStream err) throws IOException {
          mUsageTracker = new UsageTracking(new Properties(), moduleName(), CommandLine.getRunId(), configFile, true);
          return super.mainExec(out, err);
        }
      };
      final MemoryPrintStream ps = new MemoryPrintStream();
      final int[] code = new int[1];
      code[0] = 1001;
      final Runnable run = new Runnable() {
        @Override
        public void run() {
          int returncode = cli.mainInit(new String[] {"-p", "3283"}, ps.outputStream(), ps.printStream());
          synchronized (cli.mSync) {
            code[0] = returncode;
          }
        }
      };
      final Thread serverThread = new Thread(run);
      serverThread.start();
      synchronized (cli.mSync) {
        while (!cli.getStarted() && code[0] == 1001) {
          cli.mSync.wait(50);
        }
      }

      UUID runId = UUID.randomUUID();
      if (cli.getStarted()) {
        final HttpUsageTrackingClient http = new HttpUsageTrackingClient("http://localhost:" + 3283 + "/usage", new UsageConfiguration(new Properties()), false);
        runId = UUID.randomUUID();
        http.recordBeginning("testModule", runId);
        http.recordEnd(65, "testModule", runId, true);
        serverThread.interrupt();
      }
      serverThread.join();
      assertEquals(ps.toString(), 0, code[0]);
      final File[] files = dir.listFiles();
      assertNotNull(files);
      assertEquals(1, files.length);
      final File usageOut = files[0];
      final SimpleDateFormat df = new SimpleDateFormat("yyyy-MM");
      assertEquals(df.format(new Date()) + ".usage", usageOut.getName());
      final String str = FileUtils.fileToString(usageOut);
      TestUtils.containsAll(str, License.getSerialNumber(), Environment.getVersion(), "testModule", "Start", "Success", runId.toString());
      final String[] outputLines = str.split("\r\n|\n");
      final String[] line1 = outputLines[0].split("\tS="); //splits into message and md5 sum
      assertEquals(MD5Utils.md5(line1[0]), line1[1]);
      final String[] line2 = outputLines[1].split("\tS=");
      assertEquals(MD5Utils.md5(line1[1] + line2[0]), line2[1]);
    }
  }
}
