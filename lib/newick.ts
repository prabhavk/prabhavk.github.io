// lib/newick.ts
export type NwNode = {
  name?: string;
  /** Edge length from parent to this node (root has undefined) */
  length?: number;
  children?: NwNode[];
};

/** Ensure Newick ends with a single trailing semicolon. */
export function ensureSemicolon(s: string): string {
  const t = (s ?? "").trim();
  if (!t) return "";
  return t.endsWith(";") ? t : `${t};`;
}

/**
 * Minimal, robust Newick parser (names: unquoted tokens w/o ,():; or whitespace)
 * Supports edge lengths after ":"; ignores comments/whitespace.
 */
export function parseNewick(input: string): NwNode {
  if (typeof input !== "string") throw new Error("Newick must be a string");
  const s = ensureSemicolon(input);
  let i = 0;

  const isWS = (c: string) => /\s/.test(c);
  function skipWs() {
    while (i < s.length && isWS(s[i]!)) i++;
  }
  function readName(): string {
    skipWs();
    const start = i;
    while (i < s.length && !/[,\(\):;\s]/.test(s[i]!)) i++;
    return s.slice(start, i);
  }
  function readNumber(): number | undefined {
    skipWs();
    const start = i;
    while (i < s.length && /[0-9eE+\-\.]/.test(s[i]!)) i++;
    const tok = s.slice(start, i);
    if (!tok) return undefined;
    const v = Number(tok);
    return Number.isFinite(v) ? v : undefined;
  }

  function parseSubtree(): NwNode {
    skipWs();
    const node: NwNode = {};

    if (s[i] === "(") {
      i++; // consume '('
      node.children = [];
      for (;;) {
        node.children.push(parseSubtree());
        skipWs();
        if (s[i] === ",") {
          i++; // next child
          continue;
        }
        if (s[i] === ")") {
          i++; // end children
          break;
        }
        if (i >= s.length) throw new Error("Unclosed '(' in Newick");
        // tolerate stray whitespace
      }
      skipWs();
      const nm = readName();
      if (nm) node.name = nm;
    } else {
      // leaf
      const nm = readName();
      if (nm) node.name = nm;
    }

    skipWs();
    if (s[i] === ":") {
      i++;
      const len = readNumber();
      if (len !== undefined) node.length = len;
    }
    return node;
  }

  const root = parseSubtree();
  skipWs();
  if (s[i] === ";") i++;
  skipWs();
  if (i !== s.length) {
    // soft tolerance: ignore trailing junk after ';'
  }
  return root;
}
