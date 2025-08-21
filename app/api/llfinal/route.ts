// app/api/llfinal/route.ts
import { NextResponse } from "next/server";
import mysql, { RowDataPacket } from "mysql2/promise";

// Force Node runtime (mysql2 doesn't work on Edge)
export const runtime = "nodejs";
// Always return fresh DB reads
export const revalidate = 0;

type Method = "Parsimony" | "Dirichlet" | "SSH";

interface EmtrRow extends RowDataPacket {
  init_method: string;
  ll_final: number;
}

function normalizeMethod(x: unknown): Method | null {
  const s = String(x ?? "").trim().toUpperCase();
  if (s === "PARSIMONY") return "Parsimony";
  if (s === "DIRICHLET") return "Dirichlet";
  if (s === "SSH") return "SSH";
  return null;
}

/** Build a mysql2 pool explicitly instead of using a URI with ssl=… */
function makePool() {
  const { PS_HOST, PS_USERNAME, PS_PASSWORD, PS_DATABASE, DATABASE_URL } = process.env;

  if (PS_HOST && PS_USERNAME && PS_PASSWORD && PS_DATABASE) {
    return mysql.createPool({
      host: PS_HOST,
      user: PS_USERNAME,
      password: PS_PASSWORD,
      database: PS_DATABASE,
      port: 3306,
      ssl: { rejectUnauthorized: true },
      connectionLimit: 10,
    });
  }

  if (DATABASE_URL) {
    const url = new URL(DATABASE_URL);
    // Remove any unknown ssl query param; mysql2 expects ssl in options.
    url.searchParams.delete("ssl");

    return mysql.createPool({
      host: url.hostname,
      user: decodeURIComponent(url.username),
      password: decodeURIComponent(url.password),
      database: url.pathname.replace(/^\//, ""),
      port: url.port ? Number(url.port) : 3306,
      ssl: { rejectUnauthorized: true },
      connectionLimit: 10,
    });
  }

  throw new Error("Database configuration missing. Set PS_* vars or DATABASE_URL.");
}

const pool = makePool();

export async function GET(req: Request) {
  try {
    const { searchParams } = new URL(req.url);
    const root = searchParams.get("root");
    if (!root) {
      return NextResponse.json({ error: "Missing root" }, { status: 400 });
    }

    const sql = `
      SELECT init_method, ll_final
      FROM emtr_results
      WHERE root = ? AND status = 'completed'
    `;

    const [rows] = await pool.query<EmtrRow[]>(sql, [root]);

    const data = rows
      .map((r) => {
        const method = normalizeMethod(r.init_method);
        const ll_final = Number(r.ll_final);
        return method && Number.isFinite(ll_final) ? { method, ll_final } : null;
      })
      .filter(
        (x): x is { method: Method; ll_final: number } => x !== null
      );

    return NextResponse.json({ data });
  } catch (err) {
    const message = err instanceof Error ? err.message : "Internal error";
    return NextResponse.json({ error: message }, { status: 500 });
  }
}
