// @copyright 2022-2024, Antoine Basset (CNES)
// This file is part of Linx <github.com/kabasset/Linx>
// SPDX-License-Identifier: Apache-2.0

#ifndef _LINXRUN_PROGRAMOPTIONS_H
#define _LINXRUN_PROGRAMOPTIONS_H

#include "Linx/Base/TypeUtils.h" // LINX_FORWARD

#include <boost/program_options.hpp>
#include <iostream>
#include <sstream>

namespace Linx {

/// @cond
/**
 * @brief Internal shortcut.
 */
namespace po = boost::program_options;
/// @endcond

/**
 * @brief Helper class to declare positional, named and flag options, as well as some help message.
 * 
 * There are three kinds of options:
 * - positional options: see `positional()`,
 * - named options: see `named()`,
 * - flags: see `flag()`.
 * 
 * Positional and named options are optional if they are given a default value.
 * Arguments of positional options are parsed in order.
 * Named options and flags accept a short name, given after the long name, separated by a comma.
 * Arguments of named options are separated from the option name by a space or equal sign.
 * 
 * Arguments are parsed from the command line with `parse()`.
 * The function will print a nicely formatted one in case of failure or if explicitely requested, e.g. with option `--help`.
 * 
 * After parsing, arguments are queried with `as()`.
 * 
 * Here is an example command line with every kind of options:
 * 
 * `tree -d -L 2 --sort=size ~`
 * 
 * and a proposed implementation:
 * 
 * \code
 * 
 * int main(int argc, const char* argv[])
 * {
 *   ProgramOptions options("List contents of a directory");
 * 
 *   options.positional("dir", "The parent directory", std::string("."));
 *   options.flag("dirsonly,d", "List directories only");
 *   options.named("level,L", "Descend only level directories deep", 0);
 *   options.named("sort", "Select sort: name,version,size,mtime,ctime", std::string("name"));
 * 
 *   options.parse(argc, argv);
 * 
 *   const auto dir = options.as<std::string>("dir");
 *   const auto dirsonly = options.has("dirsonly");
 *   const auto level = options.as<int>("level");
 *   const auto sort = options.as<std::string>("sort");
 * 
 *   ...
 * 
 *   return 0;
 * }
 * \endcode
 */
class ProgramOptions {
private:

  /**
   * @brief Helper class to print help messages!
   */
  class Help {
  public:

    /**
     * @brief Constructor.
     */
    explicit Help(const std::string& description) :
        m_desc(description), m_usage(" [options]"), m_positionals(), m_nameds()
    {}

    /**
     * @brief Check whether an option has a short name.
     */
    static bool has_short_name(const std::string& name)
    {
      return name.length() > 3 && name[name.length() - 2] == ',';
    }

    /**
     * @brief Get the long name of an option.
     */
    static std::string long_name(const std::string& name)
    {
      if (has_short_name(name)) {
        return name.substr(0, name.length() - 2);
      }
      return name;
    }

    /**
     * @brief Declare a positional option.
     */
    void positional(const std::string& name, const std::string& description)
    {
      const auto argument = "<" + long_name(name) + ">";
      m_usage += " " + argument;
      m_positionals.emplace_back(argument + "\n      " + append_dot(description));
    }

    /**
     * @brief Declare a positional option with default value.
     */
    template <typename T>
    void positional(const std::string& name, const std::string& description, T&& default_value)
    {
      const auto argument = "<" + long_name(name) + ">";
      m_usage += " [" + argument + "]";
      m_positionals.emplace_back(argument + "\n      " + append_dot(description));
      with_default(m_positionals.back(), LINX_FORWARD(default_value));
    }

    /**
     * @brief Declare a named option.
     */
    void named(const std::string& name, const std::string& description)
    {
      auto option = has_short_name(name) ? std::string {'-', name.back(), ',', ' '} : std::string();
      const auto ln = long_name(name);
      option += "--" + ln + " <" + ln + ">\n      " + append_dot(description);
      m_nameds.push_back(std::move(option));
    }

    /**
     * @brief Declare a named option with default value.
     */
    template <typename T>
    void named(const std::string& name, const std::string& description, T&& default_value)
    {
      named(name, description);
      with_default(m_nameds.back(), LINX_FORWARD(default_value));
    }

    /**
     * @brief Declare an implicit option.
     */
    template <typename T>
    void implicit(const std::string& name, const std::string& description, T&& default_value, T&& implicit_value)
    {
      named(name, description);
      with_implicit(m_nameds.back(), LINX_FORWARD(default_value), LINX_FORWARD(implicit_value));
    }

    void flag(const std::string& name, const std::string& description)
    {
      auto option = has_short_name(name) ? std::string {'-', name.back(), ',', ' '} : std::string();
      const auto ln = long_name(name);
      option += "--" + ln + "\n      " + append_dot(description);
      m_nameds.push_back(std::move(option));
    }

    /**
     * @brief Print the help message to a given stream.
     */
    void to_stream(const std::string& argv0, std::ostream& out = std::cout)
    {
      // Help
      if (not m_desc.empty()) {
        out << "\n" << m_desc << "\n";
      }

      // Usage
      out << "\nUsage:\n\n  " << argv0 << m_usage << "\n";
      // FIXME split program name?

      // Positional options
      for (const auto& o : m_positionals) {
        out << "\n  " << o;
      }
      if (not m_positionals.empty()) {
        out << "\n";
      }

      // Named options
      if (m_nameds.empty()) {
        return;
      }
      out << "\nOptions:\n";
      for (const auto& o : m_nameds) {
        out << "\n  " << o;
      }

      out << "\n\n";
      std::flush(out);
    }

  private:

    /**
     * @brief Add a default value to an option.
     */
    template <typename T>
    void with_default(std::string& option, T&& value)
    {
      if constexpr (std::is_same_v<std::decay_t<T>, char>) {
        option.append("\n      [default: " + std::string {value} + "]");
      } else if constexpr (std::is_same_v<std::decay_t<T>, std::string>) {
        option.append("\n      [default: " + LINX_FORWARD(value) + "]");
      } else {
        option.append("\n      [default: " + std::to_string(LINX_FORWARD(value)) + "]");
      }
    }

    /**
     * @brief Add default and implicit values to an option.
     */
    template <typename T>
    void with_implicit(std::string& option, T&& default_value, T&& implicit_value)
    {
      with_default(option, default_value);
      if constexpr (std::is_same_v<std::decay_t<T>, char>) {
        option.append("\n      [implicit: " + std::string {implicit_value} + "]");
      } else if constexpr (std::is_same_v<std::decay_t<T>, std::string>) {
        option.append("\n      [implicit: " + LINX_FORWARD(implicit_value) + "]");
      } else {
        option.append("\n      [implicit: " + std::to_string(LINX_FORWARD(implicit_value)) + "]");
      }
    }

    /**
     * @brief Append a dot if there is not.
     */
    static std::string append_dot(const std::string description)
    {
      if (description.back() == '.') { // FIXME other punctuation marks
        return description;
      }
      return description + '.';
    }

  private:

    std::string m_desc; ///< The program description
    std::string m_usage; ///< The program usage
    std::vector<std::string> m_positionals; ///< The positional options description
    std::vector<std::string> m_nameds; ///< The named options description
  };

public:

  /**
   * @brief Make a `ProgramOptions` with optional description string and help option.
   * @param description The program description
   * @param help The help option (disable with empty string)
   */
  ProgramOptions(const std::string& description = "", const std::string& help = "help,h") :
      m_named("Options", 120), m_add(m_named.add_options()), m_positional(), m_variables(), m_desc(description),
      m_help(help)
  {
    if (m_help.length() > 0) {
      flag(m_help.c_str(), "Print help message");
      m_help = Help::long_name(m_help);
    }
  }

  /**
   * @brief Declare a positional option.
   */
  template <typename T>
  void positional(const char* name, const char* description)
  {
    positional(name, po::value<T>()->required(), description);
    m_desc.positional(name, description);
  }

  /**
   * @brief Declare a positional option with default value.
   */
  template <typename T>
  void positional(const char* name, const char* description, T default_value)
  {
    positional(name, po::value<T>()->default_value(default_value), description);
    m_desc.positional(name, description, default_value);
  }

  /**
   * @brief Declare a named option.
   */
  template <typename T>
  void named(const char* name, const char* description)
  {
    named(name, po::value<T>()->required(), description);
    m_desc.named(name, description);
  }

  /**
   * @brief Declare a named option with default value.
   */
  template <typename T>
  void named(const char* name, const char* description, T default_value)
  {
    named(name, po::value<T>()->default_value(default_value), description);
    m_desc.named(name, description, default_value);
  }

  /**
   * @brief Declare an implicit option.
   * 
   * An implicit option has two default values:
   * the first one is used when the option is not there, and
   * the second one is used when the option is there with no argument.
   */
  template <typename T>
  void implicit(const std::string& name, const std::string& description, T&& default_value, T&& implicit_value)
  {
    named(name, po::value<T>()->default_value(default_value)->implicit_value(implicit_value), description);
    m_desc.implicit(name, description, default_value, implicit_value);
  }

  /**
   * @brief Declare a flag option.
   */
  void flag(const char* name, const char* description)
  {
    named(name, po::value<bool>()->default_value(false)->implicit_value(true), description);
    m_desc.flag(name, description);
  }

  /**
   * @brief Parse a command line.
   * 
   * If the help option was enabled and is in the command line, then the help message is printed and the program completes.
   * If the parsing fails, then the help message is printed to the standard error and the program terminates with error code.
   */
  void parse(int argc, const char* const argv[])
  {
    try {
      po::store(po::command_line_parser(argc, argv).options(m_named).positional(m_positional).run(), m_variables);
      if (not m_help.empty() && has(m_help.c_str())) {
        m_desc.to_stream(argv[0]);
        exit(0);
      }
      po::notify(m_variables);
    } catch (...) {
      std::cerr << "\nFATAL: Cannot parse command line.\n";
      m_desc.to_stream(argv[0], std::cerr);
      std::rethrow_exception(std::current_exception());
    }
  }

  /**
   * @brief Parse a command line (vector of `const char*` arguments).
   */
  void parse(const std::vector<const char*>& args)
  {
    parse(args.size(), args.data());
  }

  /**
   * @brief Parse a command line (vector of string arguments).
   */
  void parse(const std::vector<std::string>& args)
  {
    std::vector<const char*> cstr(args.size());
    std::transform(args.begin(), args.end(), cstr.begin(), [](const auto& a) {
      return a.c_str();
    });
    parse(cstr);
  }

  /**
   * @brief Parse a command line (space-separated arguments as a single string).
   */
  void parse(const std::string& args)
  {
    std::istringstream iss(args);
    parse(std::vector<std::string>(std::istream_iterator<std::string> {iss}, std::istream_iterator<std::string>()));
  }

  /**
   * @brief Check whether a given option is set.
   */
  bool has(const char* name) const
  {
    try {
      return as<bool>(name); // Throw for non booleans
    } catch (boost::bad_any_cast&) {
      return m_variables.count(name); // Incompatible with flags
    }
  }

  /**
   * @brief Get the value of a given option.
   * 
   * Throws if the option is not set.
   */
  template <typename T>
  T as(const char* name) const
  {
    return m_variables[name].as<T>();
  }

private:

  /**
   * @brief Declare a positional option with custom semantics.
   */
  void positional(const char* name, const po::value_semantic* value, const char* description)
  {
    m_add(name, value, description);
    const int max_args = value->max_tokens();
    m_positional.add(name, max_args);
  }

  /**
   * @brief Declare a named option with custom semantics.
   * 
   * A short form (1-character) of the option can be provided, separated by a comma.
   */
  void named(const char* name, const po::value_semantic* value, const char* description)
  {
    m_add(name, value, description);
  }

private:

  po::options_description m_named;
  po::options_description_easy_init m_add;
  po::positional_options_description m_positional;
  po::variables_map m_variables;
  Help m_desc;
  std::string m_help;
};

} // namespace Linx

#endif
