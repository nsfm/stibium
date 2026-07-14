#ifndef VM_IMPORT_H
#define VM_IMPORT_H

#include <QString>

/*
 *  Implements --import-vm: translates a Fidget .vm math tape
 *  (github.com/mkeeter/fidget) into a Stibium .sb project holding a
 *  single shape node.  See doc/FOREIGN-IMPORT.md for the format
 *  analysis this follows.
 *
 *  Returns a process exit code (0 on success); prints diagnostics to
 *  stderr.
 */
int importVmHeadless(const QString& vm_path, const QString& sb_path);

#endif
