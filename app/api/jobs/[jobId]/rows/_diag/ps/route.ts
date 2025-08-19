import { NextResponse } from "next/server";
import { connect } from "@planetscale/database";

export const runtime = "edge";

export async function GET() {
  try {
    const host = process.env.PS_HOST;
    const username = process.env.PS_USERNAME;
    const password = process.env.PS_PASSWORD;

    if (!host || !username || !password) {
      return NextResponse.json(
        {
          ok: false,
          step: "env",
          hasHost: !!host,
          hasUsername: !!username,
          hasPassword: !!password,
        },
        { status: 500 }
      );
    }

    const conn = connect({ host, username, password });
    const r = await conn.execute("SELECT 1 AS one");
    return NextResponse.json({ ok: true, result: r.rows });
  } catch (e) {
    return NextResponse.json(
      {
        ok: false,
        step: "query",
        error: String(e),
      },
      { status: 500 }
    );
  }
}
