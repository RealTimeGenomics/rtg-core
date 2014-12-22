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

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

import com.rtg.util.MD5Utils;
import com.rtg.util.StringUtils;

/**
 * container for usage message
 * <p>lines look like</p>
 * <p><code>timestamp\tlicense_serial\trun_uuid\trtg_version\tmodule_name\trun_status\trun_metric\tuser_name\thost_name\tcommand_line\tchecksum</code></p>
 * <p><code>run_metric</code> only contains a value on end messages</p>. <code>run_status</code> typically should be Start/Success/Fail. <code>user_name/host_name/command_line</code> are optional based on the users configuration and need not be set, in this case they will take the value <code>N/A</code>.
 */
public final class UsageMessage {
  static final String SIGNATURE = "S=";
  /** length allowed for user name in database */
  public static final int USERNAME_TRIM_LENGTH = 30;
  /** silly length allowed for host name in database */
  public static final int HOSTNAME_TRIM_LENGTH = 30;
  /** length allowed for command line in database */
  public static final int COMMANDLINE_TRIM_LENGTH = 1000;

  private String mSerialNumber;
  private String mRunId;
  private String mVersion;
  private String mModule;
  private String mMessageType;
  private String mMetric = "N/A";
  private String mDate;

  String mUsername = "N/A";
  String mHostname = "N/A";
  private String mCommandLine = "N/A";

  private String mChecksum;

  private UsageMessage() { }

  /**
   * helper for <code>start</code> message
   * @param serialNumber license serial number
   * @param runId unique id for run
   * @param version version string
   * @param module module being run
   * @return the usage message container
   */
  static UsageMessage startMessage(String serialNumber, String runId, String version, String module) {
    final UsageMessage ret = new UsageMessage();
    ret.setSerialNumber(serialNumber);
    ret.setRunId(runId);
    ret.setVersion(version);
    ret.setModule(module);
    ret.setMessageType("Start");
    return ret;
  }

  /**
   * helper for message
   * @param serialNumber license serial number
   * @param runId unique id for run
   * @param version version string
   * @param module module being run
   * @param metric usage metric for run
   * @param status should be &quot;Success&quot; for a successful run and &quot;Fail&quot; for a failed run
   * @return the message
   */
  static UsageMessage setMessage(String serialNumber, String runId, String version, String module, String metric, String status) {
    final UsageMessage ret = new UsageMessage();
    ret.setSerialNumber(serialNumber);
    ret.setRunId(runId);
    ret.setVersion(version);
    ret.setModule(module);
    ret.setMessageType(status);
    if (metric != null) {
      ret.setMetric(metric);
    }
    return ret;
  }

  public void setSerialNumber(String serialNumber) {
    mSerialNumber = serialNumber;
  }

  public void setRunId(String runId) {
    mRunId = runId;
  }

  public void setVersion(String version) {
    mVersion = version;
  }

  public void setModule(String module) {
    mModule = module;
  }

  public void setMessageType(String messageType) {
    mMessageType = messageType;
  }

  public void setMetric(String metric) {
    mMetric = metric;
  }

  /**
   * sets the date using date format <code>yyyy-MM-dd HH:mm:ss</code>
   * @param date value to set date to
   */
  public void setDate(Date date) {
    final DateFormat messageDateFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
    mDate = messageDateFormat.format(date);
  }

  /**
   * @param username value to set, will truncate if needed
   */
  public void setUsername(String username) {
    if (username != null) {
      mUsername = trimField(username, USERNAME_TRIM_LENGTH);
    }
  }

  /**
   * @param hostname value to set, will truncate if needed
   */
  public void setHostname(String hostname) {
    if (hostname != null) {
      mHostname = trimField(hostname, HOSTNAME_TRIM_LENGTH);
    }
  }

  /**
   * sets the command line, if value is longer than 1000 characters it is truncated at that length.
   * replaces any whitespace with space characters
   * @param commandLine value to set command line to
   */
  public void setCommandLine(String commandLine) {
    if (commandLine != null) {
      mCommandLine = trimField(commandLine, COMMANDLINE_TRIM_LENGTH);
    }
  }

  /**
   * if value is longer than {@code trimLength} characters it is truncated at that length.
   * replaces any whitespace with space characters
   * @param fieldValue value to set command line to
   * @param trimLength length to truncate to
   * @return truncated and space replaced string
   */
  public static String trimField(String fieldValue, int trimLength) {
    if (fieldValue == null) {
      return null;
    }
    String ret;
    if (fieldValue.length() > trimLength) {
      ret = fieldValue.substring(0, trimLength);
    } else {
      ret = fieldValue;
    }
    ret = ret.replaceAll("\\s", " ");
    return ret;
  }

  /**
   * get the checksum for this message, only available after {@link UsageMessage#formatLine(String)} has been called.
   * @return the checksum
   */
  public String getChecksum() {
    return mChecksum;
  }

  /**
   * formats the message for writing to file and sets the value for the {@link com.rtg.usage.UsageMessage#getChecksum()} method
   * @param prevChecksum checksum for last message, or null if this is the first message
   * @return the formatted message (including new line)
   */
  public String formatLine(String prevChecksum) {
    final StringBuilder messageBuilder = new StringBuilder();
    messageBuilder.append(mDate); //Date first
    final String[] fields = {mSerialNumber, mRunId, mVersion, mModule, mMessageType, mMetric, mUsername, mHostname, mCommandLine}; //rest in this order
    for (String s : fields) {
      messageBuilder.append("\t").append(s.replaceAll("\\s", " "));
    }
    final String messageContents = messageBuilder.toString();
    mChecksum = sign((prevChecksum == null ? "" : prevChecksum) + messageContents);
    messageBuilder.append("\t").append(SIGNATURE).append(mChecksum);
    messageBuilder.append(StringUtils.LS);
    return messageBuilder.toString();
  }

  static String sign(String message) {
    return MD5Utils.md5(message);
  }
}
