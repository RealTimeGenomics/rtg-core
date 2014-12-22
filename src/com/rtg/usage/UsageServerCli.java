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
import java.net.URI;
import java.net.URISyntaxException;

import com.rtg.launcher.AbstractCli;
import com.rtg.util.StringUtils;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

/**
 * module entrance for usage server
 */
public class UsageServerCli extends AbstractCli {

  private static final String PORT = "port";
  private static final String THREADS = "threads";

  // these can be used to check if the server has been started if we're called programmatically
  protected final Object mSync;
  private boolean mStarted = false;

  /**
   * Constructor
   */
  public UsageServerCli() {
    mSync = new Object();
    mSuppressUsage = true;
  }

  boolean getStarted() {
    synchronized (mSync) {
      return mStarted;
    }
  }

  @Override
  protected void initFlags() {
    final String host = mUsageTracker.getUsageConfiguration().getUsageHost();
    int defaultPort = 8080;
    if (host != null) {
      try {
        final URI uri = new URI(host);
        if (uri.getScheme() != null && uri.getScheme().equalsIgnoreCase("http")) {
          final int hostPort = uri.getPort();
          defaultPort = hostPort == -1 ? 80 : hostPort;
        }
      } catch (URISyntaxException e) {
        throw new NoTalkbackSlimException("Malformed usage host specified in usage configuration: " + host);
      }
    }
    CommonFlagCategories.setCategories(mFlags);
    mFlags.setDescription("Start usage tracking server.");
    mFlags.registerOptional('p', PORT, Integer.class, "INT", "port on which to listen for usage logging connections.", defaultPort).setCategory(CommonFlagCategories.UTILITY);
    mFlags.registerOptional('T', THREADS, Integer.class, "INT", "number of worker threads to handle incoming connections.", 4).setCategory(CommonFlagCategories.UTILITY);
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    Diagnostic.info("Checking usage configuration.");
    if (!mUsageTracker.getUsageConfiguration().isEnabled()) {
      throw new NoTalkbackSlimException("Cannot start usage server without RTG_USAGE configuration option set. " + UsageTracking.SEE_MANUAL);
    } else if (mUsageTracker.getUsageConfiguration().getUsageDir() == null) {
      throw new NoTalkbackSlimException("Cannot start usage server without RTG_USAGE_DIR configuration option set. " + UsageTracking.SEE_MANUAL);
    }
    final Integer port = (Integer) mFlags.getValue(PORT);
    final String configHost = mUsageTracker.getUsageConfiguration().getUsageHost();
    // Output some warnings if the port doesn't correspond with where the clients will be pointing
    if (configHost != null) {
      try {
        final URI uri = new URI(configHost);
        if (uri.getScheme() != null && uri.getScheme().equalsIgnoreCase("http")) {
          final int configPort = uri.getPort() == -1 ? 80 : uri.getPort();
          if (configPort != port) {
            Diagnostic.warning("Specified port " + port + " does not correspond with port from usage configuration: " + configPort);
          }
        } else {
          Diagnostic.warning("This server (HTTP on port " + port + ") does not correspond to current usage configuration: " + configHost);
        }
      } catch (URISyntaxException e) {
        throw new NoTalkbackSlimException("Malformed usage host URI specified in usage configuration: " + configHost);
      }
    } else {
      Diagnostic.warning("Clients will not be able to connect until RTG_USAGE_HOST has been set. " + UsageTracking.SEE_MANUAL);
    }
    final UsageServer us = new UsageServer(port, new File(mUsageTracker.getUsageConfiguration().getUsageDir()), (Integer) mFlags.getValue(THREADS));
    synchronized (mSync) {
      us.start();
      out.write(("Usage server listening on port " + us.getPort() + StringUtils.LS).getBytes());
      mStarted = true;
      try {
        boolean cont = true;
        while (cont) {
          try {
            mSync.wait(1000);
          } catch (InterruptedException e) {
            cont = false;
          }
        }
      } finally {
        us.end();
      }
    }
    return 0;
  }

  @Override
  public String moduleName() {
    return "usageserver";
  }
}
