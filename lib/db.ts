// lib/db.ts
import mysql from "mysql2/promise";

export const pool = mysql.createPool({
  host: process.env.EMTR_DB_HOST ?? "127.0.0.1",
  user: process.env.EMTR_DB_USER ?? "emtr",
  password: process.env.EMTR_DB_PASS ?? "emtrpw",
  database: process.env.EMTR_DB_NAME ?? "emtr",
  port: Number(process.env.EMTR_DB_PORT ?? 3306),
  waitForConnections: true,
  connectionLimit: 10,
  maxIdle: 10,
  idleTimeout: 60000,
  queueLimit: 0,
});
