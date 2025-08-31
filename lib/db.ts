// lib/db.ts
import mysql from "mysql2/promise";

export type SQLParam = string | number | boolean | Date | null | undefined;
type Row = Record<string, unknown>;

/** Tell TypeScript about our global pool slot (no eslint-disable needed). */
declare global {
  // This is a declaration merge on the global scope; it does not emit JS.
  // eslint-disable-next-line no-var
  var __MYSQL_POOL__: mysql.Pool | undefined;
}

/** Create a new pool (PlanetScale requires TLS). We ignore any ?ssl=... in the URL and set SSL in code. */
function createPool(): mysql.Pool {
  const url = process.env.DATABASE_URL;
  if (!url) throw new Error("DATABASE_URL is not set");

  const u = new URL(url);
  const host = u.hostname;
  const port = u.port ? Number(u.port) : 3306;
  const user = decodeURIComponent(u.username);
  const password = decodeURIComponent(u.password);
  const database = u.pathname.replace(/^\//, "");

  return mysql.createPool({
    host,
    port,
    user,
    password,
    database,
    ssl: { rejectUnauthorized: true }, // PlanetScale/Vitess-friendly
    waitForConnections: true,
    connectionLimit: 10,
  });
}

/** Reuse a single pool across HMR in dev. */
function getPool(): mysql.Pool {
  if (!globalThis.__MYSQL_POOL__) {
    globalThis.__MYSQL_POOL__ = createPool();
  }
  return globalThis.__MYSQL_POOL__;
}

/** Run a parameterized SELECT and return typed rows. */
export async function q<T extends Row = Row>(
  sql: string,
  params: ReadonlyArray<SQLParam> = [],
): Promise<T[]> {
  const [rows] = await getPool().execute<mysql.RowDataPacket[]>(
    sql,
    params as unknown[],
  );
  return rows as unknown as T[];
}

/** Like q(), but returns the first row (or null). */
export async function q1<T extends Row = Row>(
  sql: string,
  params: ReadonlyArray<SQLParam> = [],
): Promise<T | null> {
  const rows = await q<T>(sql, params);
  return rows.length ? rows[0] : null;
}

/** Execute a non-SELECT (INSERT/UPDATE/DELETE) and return the result header. */
export async function exec(
  sql: string,
  params: ReadonlyArray<SQLParam> = [],
): Promise<mysql.ResultSetHeader> {
  const [res] = await getPool().execute<mysql.ResultSetHeader>(
    sql,
    params as unknown[],
  );
  return res;
}

/** Back-compat alias if other code imported `query` previously. */
export const query = q;

/** Optional: expose the pool getter if you ever need raw connections/transactions. */
export { getPool };
