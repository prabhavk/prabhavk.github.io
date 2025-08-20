// lib/db.ts
import mysql from "mysql2/promise";

export type SQLParam = string | number | boolean | Date | null | undefined;
type Row = Record<string, unknown>;

/** Tell TypeScript about our global pool slot (no eslint-disable needed). */
declare global {
  // This is a declaration merge on the global scope; it does not emit JS.
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
    ssl: { rejectUnauthorized: true },
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

/** Run a parameterized query and return rows typed as T. */
export async function query<T extends Row = Row>(
  sql: string,
  params: ReadonlyArray<SQLParam> = [],
): Promise<T[]> {
  const [rows] = await getPool().execute<mysql.RowDataPacket[]>(
    sql,
    params as unknown[],
  );
  return rows as unknown as T[];
}
