declare module '@bunchmark/stats' {
  export type MWUResult = { U: number; z: number; p: number };
  export function mwu(a: readonly number[], b: readonly number[]): MWUResult;
}
  