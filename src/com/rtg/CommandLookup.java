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

import java.util.Locale;
import java.util.NavigableSet;
import java.util.TreeSet;

import com.rtg.util.StringUtils;

/**
 * Allows lookup of all commands that are part of a product and retrieval of a command by name.
 */
public abstract class CommandLookup {

  /**
   * @return list of the enum values for commands
   */
  abstract Command[] commands();

  /**
   * Look a command by name, without prefix expansion to licensed commands.
   * @param name the input command name or prefix
   * @return the matching command, or null if an invalid command name was given
   */
  public Command findModule(final String name) {
    final String norm = name.toUpperCase(Locale.getDefault());
    for (Command module : commands()) {
      if (module.getCommandName().equals(norm)) {
        return module;
      }
    }
    return null;
  }

  /**
   * Look a command by name, allowing prefix expansion to licensed commands.
   * @param name the input command name or prefix
   * @return the matching command, or null if no command name matched
   */
  public Command findModuleWithExpansion(final String name) {
    final String norm = name.toUpperCase(Locale.getDefault());

    final NavigableSet<String> names = new TreeSet<>();
    for (Command module : commands()) {
      if (module.isLicensed()) {
        names.add(module.getCommandName());
      }
    }

    return findModule(StringUtils.expandPrefix(names, norm));
  }
}
