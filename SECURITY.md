# Security Policy

## Supported Versions

We maintain security updates for the following versions:

| Version | Supported          |
|---------|--------------------|
| Current | ✅ Yes             |
| < 1.0   | ✅ Yes             |

## Reporting a Vulnerability

### Data Privacy & Local Execution

**Important**: This tool runs entirely on your local machine and does NOT upload any data to external servers. All processing happens locally within your environment.

### How to Report

If you discover a security vulnerability, please report it privately rather than creating a public issue.

**Steps to report:**

1. **Email**: Send details to [myzhang0726@foxmail.com]
2. **Format**: Include "SECURITY" in the subject line
3. **Information to include**:
   - Description of the vulnerability
   - Steps to reproduce
   - Potential impact
   - Suggested fix (if known)

### What We Consider a Vulnerability

We consider the following as security issues that should be reported:

- **Local data exposure**: Unintended access to files on your system
- **Command injection**: Ability to execute arbitrary commands
- **Path traversal**: Access to files outside intended directories
- **Local file corruption**: Unintended modification of local files
- **Information leakage**: Sensitive data exposed in logs or error messages
- **Denial of service**: Tool crashes or hangs due to malformed input

### What We DON'T Consider a Vulnerability

Given this is a local-only tool:

- Access to files in directories you explicitly provide
- Operations that require explicit user permission
- Resource usage that matches user expectations
- Issues requiring physical access to your machine

### Response Timeline

- **Initial acknowledgment**: Within 48 hours
- **Investigation**: 3-7 business days
- **Fix & release**: As soon as practical, typically within 30 days

### After Reporting

- We will confirm receipt of your report
- We will keep you informed of our progress
- We will credit you in the release notes (if desired)
- We will maintain confidentiality until the fix is released

## Best Practices for Users

1. **Verify downloads**: Only download from official releases
2. **Check permissions**: Review what files the tool accesses
3. **Keep updated**: Use the latest version for security patches
4. **Backup data**: Always backup important data before using new tools
5. **Review source**: Check the code for peace of mind

## Security Principles

- **Local-first**: No data leaves your machine
- **Minimal permissions**: Only accesses what's necessary
- **Transparent**: Source code is available for review
- **No telemetry**: No usage data is collected

## Contact

For security-related questions not related to vulnerability reports, please open a discussion or issue on GitHub.

---

*Last updated: April 2026*
