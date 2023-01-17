# from sqlalchemy import *
# from sqlalchemy import event
# from sqlalchemy.engine import Engine
# from sqlalchemy.orm import relationship, declarative_base, Session
#
#
# def connect_network_db(db_file: str) -> Engine:
#     return create_engine("sqlite+pysqlite:///"+db_file, echo=False, future=True)
#
#
# Base = declarative_base()
#
# #
# # Commented columns are columns from PyPSA that are not used
# #
#
#
# class Bus(Base):
#     __tablename__ = 'buses'
#
#     id = Column("Bus", Text, primary_key=True, index=True)
#     v_nom = Column(Float)
#     symbol = Column(Text)
#     under_construction = Column(Boolean)
#     # tags = Column(Text)
#     x = Column(Float)
#     y = Column(Float)
#     carrier = Column(Text)
#     country = Column(Text)
#     # substation_lv = Column(Boolean)
#     # generator_substation_lv = Column(Boolean)  # updated PyPSA-Eur code to extend substation_lv
#     # substation_off = Column(Boolean)
#     # type = Column(Text)
#     # unit = Column(Text)
#     # v_mag_pu_set = Column(Float)
#     # v_mag_pu_min = Column(Float)
#     # v_mag_pu_max = Column(Float)
#     # control = Column(Text)
#     # sub_network = Column(Text)
#     in_synchronous_network = Column(Boolean)
#
#     generators = relationship("Generator", back_populates="bus")
#     storage_units = relationship("StorageUnit", back_populates="bus")
#     loads = relationship("Load", back_populates="bus")
#     shunt_impedances = relationship("ShuntImpedance", back_populates="bus")
#     external_links = relationship("ExternalLink", back_populates="bus")
#
#     # only define 1 connection side of branches
#     lines = relationship("Line", back_populates="bus0", foreign_keys='Line.bus0_id')
#     transformers = relationship("Transformer", back_populates="bus0", foreign_keys='Transformer.bus0_id')
#     links = relationship("Link", back_populates="bus0", foreign_keys='Link.bus0_id')
#
#
# class Generator(Base):
#     __tablename__ = 'generators'
#
#     id = Column('Generator', Text, primary_key=True, index=True)
#     name = Column(Text)
#     carrier = Column(Text)
#     bus_id = Column('bus', ForeignKey('buses.Bus'))
#     # p_nom_min = Column(Float)
#     p_nom = Column(Float)
#     # p_nom_extendable = Column(Boolean)
#     # efficiency = Column(Float)
#     # marginal_cost = Column(Float)
#     # capital_cost = Column(Float)
#     # build_year = Column(Integer)
#     # lifetime = Column(Float)
#     # p_nom_max = Column(Float)
#     p_max_pu = Column(Float)
#     # min_up_time = Column(Integer)
#     # min_down_time = Column(Integer)
#     # up_time_before = Column(Integer)
#     # down_time_before = Column(Integer)
#     # weight = Column(Float)
#     # control = Column(Text)
#     # type = Column(Text)
#     p_min_pu = Column(Float)
#     p_set = Column(Float)
#     q_set = Column(Float)
#     # sign = Column(Float)
#     # committable = Column(Boolean)
#     # start_up_cost = Column(Float)
#     # shut_down_cost = Column(Float)
#     # ramp_limit_up = Column(Float)
#     # ramp_limit_down = Column(Float)
#     # ramp_limit_start_up = Column(Float)
#     # ramp_limit_shut_down = Column(Float)
#     # p_nom_opt = Column(Float)
#
#     bus = relationship('Bus', back_populates="generators")
#
#
# class Load(Base):
#     __tablename__ = 'loads'
#
#     id = Column('Load', Text, primary_key=True, index=True)
#     bus_id = Column('bus', ForeignKey('buses.Bus'))
#     # carrier = Column(Text)
#     # type = Column(Text)
#     p_set = Column(Float)
#     q_set = Column(Float)
#     # sign = Column(Float)
#
#     bus = relationship('Bus', back_populates="loads")
#
#
# class ShuntImpedance(Base):
#     __tablename__ = 'shunt_impedances'
#
#     id = Column('ShuntImpedance', Text, primary_key=True, index=True)
#     bus_id = Column('bus', ForeignKey('buses.Bus'))
#     # g = Column(Float)
#     # b = Column(Float)
#     # sign = Column(Float)
#     g_pu = Column(Float)
#     b_pu = Column(Float)
#
#     bus = relationship('Bus', back_populates="shunt_impedances")
#
#
# class StorageUnit(Base):
#     __tablename__ = 'storage_units'
#
#     id = Column('StorageUnit', Text, primary_key=True, index=True)
#     carrier = Column(Text)
#     bus_id = Column('bus', ForeignKey('buses.Bus'))
#     p_nom = Column(Float)
#     # capital_cost = Column(Float)
#     # max_hours = Column(Float)
#     # efficiency_store = Column(Float)
#     # efficiency_dispatch = Column(Float)
#     # cyclic_state_of_charge = Column(Boolean)
#     p_min_pu = Column(Float)
#     p_max_pu = Column(Float)
#     # build_year = Column(Integer)
#     # control = Column(Text)
#     # type = Column(Text)
#     # p_nom_extendable = Column(Boolean)
#     # p_nom_min = Column(Float)
#     # p_nom_max = Column(Float)
#     p_set = Column(Float)
#     q_set = Column(Float)
#     # sign = Column(Float)
#     # marginal_cost = Column(Float)
#     # lifetime = Column(Float)
#     # state_of_charge_initial = Column(Float)
#     # state_of_charge_initial_per_period = Column(Boolean)
#     # state_of_charge_set = Column(Float)
#     # cyclic_state_of_charge_per_period = Column(Boolean)
#     # standing_loss = Column(Float)
#     # inflow = Column(Float)
#     # p_nom_opt = Column(Float)
#
#     bus = relationship('Bus', back_populates="storage_units")
#
#
# class Line(Base):
#     __tablename__ = 'lines'
#
#     id = Column('Line', Text, primary_key=True, index=True)
#     bus0_id = Column("bus0", ForeignKey('buses.Bus'))
#     bus1_id = Column("bus1", ForeignKey('buses.Bus'))
#     # v_nom = Column(BigInteger)
#     num_parallel = Column(Float)
#     length = Column(Float)
#     # underground = Column(Boolean)
#     under_construction = Column(Boolean)
#     # tags = Column(Text)
#     # geometry = Column(Text)
#     # type = Column(Text)
#     s_max_pu = Column(Float)
#     s_nom = Column(Float)
#     # capital_cost = Column(Float)
#     # build_year = Column(Integer)
#     # x = Column(Float)
#     # r = Column(Float)
#     # g = Column(Float)
#     # b = Column(Float)
#     # s_nom_extendable = Column(Boolean)
#     # s_nom_min = Column(Float)
#     # s_nom_max = Column(Float)
#     # lifetime = Column(Float)
#     # carrier = Column(Text)
#     # terrain_factor = Column(Float)
#     # v_ang_min = Column(Float)
#     # v_ang_max = Column(Float)
#     # sub_network = Column(Text)
#     x_pu = Column(Float)
#     r_pu = Column(Float)
#     g_pu = Column(Float)
#     b_pu = Column(Float)
#     # x_pu_eff = Column(Float)
#     # r_pu_eff = Column(Float)
#     # s_nom_opt = Column(Float)
#
#     bus0 = relationship('Bus', foreign_keys=[bus0_id])
#     bus1 = relationship('Bus', foreign_keys=[bus1_id])
#
#
# class Link(Base):
#     __tablename__ = 'links'
#
#     id = Column('Link', Text, primary_key=True, index=True)
#     bus0_id = Column("bus0", ForeignKey('buses.Bus'))
#     bus1_id = Column("bus1", ForeignKey('buses.Bus'))
#     # carrier = Column(Text)
#     # geometry = Column(Text)
#     length = Column(Float)
#     p_nom = Column(Float)
#     # tags = Column(Text)
#     under_construction = Column(Boolean)
#     # underground = Column(Boolean)
#     # p_min_pu = Column(Float)
#     # underwater_fraction = Column(Float)
#     # build_year = Column(Integer)
#     # capital_cost = Column(Float)
#     # type = Column(Text)
#     # efficiency = Column(Float)
#     # lifetime = Column(Float)
#     # p_nom_extendable = Column(Boolean)
#     # p_nom_min = Column(Float)
#     # p_nom_max = Column(Float)
#     p_set = Column(Float)
#     # p_max_pu = Column(Float)
#     # marginal_cost = Column(Float)
#     # terrain_factor = Column(Float)
#     # ramp_limit_up = Column(Float)
#     # ramp_limit_down = Column(Float)
#     # p_nom_opt = Column(Float)
#
#     bus0 = relationship('Bus', foreign_keys=[bus0_id])
#     bus1 = relationship('Bus', foreign_keys=[bus1_id])
#
#
# class Transformer(Base):
#     __tablename__ = 'transformers'
#
#     id = Column('Transformer', Text, primary_key=True, index=True)
#     bus0_id = Column("bus0", ForeignKey('buses.Bus'))
#     bus1_id = Column("bus1", ForeignKey('buses.Bus'))
#     s_nom = Column(Float)
#     # tap_side = Column(Integer)
#     # tap_position = Column(Integer)
#     # build_year = Column(Integer)
#     # type = Column(Text)
#     # model = Column(Text)
#     # x = Column(Float)
#     # r = Column(Float)
#     # g = Column(Float)
#     # b = Column(Float)
#     # s_nom_extendable = Column(Boolean)
#     # s_nom_min = Column(Float)
#     # s_nom_max = Column(Float)
#     # s_max_pu = Column(Float)
#     # capital_cost = Column(Float)
#     # num_parallel = Column(Float)
#     # tap_ratio = Column(Float)
#     # phase_shift = Column(Float)
#     # lifetime = Column(Float)
#     # v_ang_min = Column(Float)
#     # v_ang_max = Column(Float)
#     # sub_network = Column(Text)
#     x_pu = Column(Float)
#     r_pu = Column(Float)
#     g_pu = Column(Float)
#     b_pu = Column(Float)
#     # x_pu_eff = Column(Float)
#     # r_pu_eff = Column(Float)
#     # s_nom_opt = Column(Float)
#
#     bus0 = relationship('Bus', foreign_keys=[bus0_id])
#     bus1 = relationship('Bus', foreign_keys=[bus1_id])
#
#
# # for cross-border flows to outside the nordics
# class ExternalLink(Base):
#     __tablename__ = 'external_links'
#
#     id = Column('ExternalLink', Integer, primary_key=True, index=True)
#     bus_id = Column('bus', ForeignKey('buses.Bus'))
#     bus_outside_id = Column('bus_outside', Text)
#     link = Column(Text)
#     p_set = Column(Float)
#     q_set = Column(Float, default=0)
#
#     bus = relationship('Bus', back_populates="external_links")
#
#
# # filter database based on column definitions in _database.py
# # it deletes columns that are not defined in _tables to get a cleaner database which only contains 'wanted' data
# # based on https://stackoverflow.com/questions/70392123/copy-one-database-to-another-using-sqlalchemy
# def filter_database(db_raw: str, db: str):
#     open(db, 'w').close()  # create empty sqlite
#
#     src_engine = connect_network_db(db_raw)
#     src_metadata = MetaData(bind=src_engine)
#
#     tgt_engine = connect_network_db(db)
#     tgt_metadata = MetaData(bind=tgt_engine)
#
#     @event.listens_for(src_metadata, "column_reflect")
#     def genericize_datatypes(inspector, tablename, column_dict):
#         column_dict["type"] = column_dict["type"].as_generic(allow_nulltype=True)
#
#     src_session = Session(src_engine)
#     tgt_session = Session(tgt_engine)
#     tgt_metadata.reflect()
#
#     tgt_metadata.clear()
#     tgt_metadata.reflect()
#     src_metadata.reflect()
#
#     # create all tables in target database
#     Base.metadata.create_all(tgt_engine)
#
#     # refresh metadata before you can copy data
#     tgt_metadata.clear()
#     tgt_metadata.reflect()
#
#     # Copy all data from src to target
#     all_tables = Base.__subclasses__()
#     tables = [table for table in all_tables if table.__tablename__ in src_metadata.tables]
#     for table in tables:
#         # data = src_conn.execute(select(src_table))
#
#         data = src_session.scalars(
#             select(table)
#         ).all()
#
#         for row in data:
#             row = row.__dict__
#             row.pop('_sa_instance_state')
#
#             new = table(**row)
#             tgt_session.add(new)
#
#         tgt_session.flush()
#     tgt_session.commit()
#     src_session.close()
#     tgt_session.close()
